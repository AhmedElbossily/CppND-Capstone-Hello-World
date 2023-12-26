#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

#include <matplot/matplot.h>

#include <cmath>

#include "util.h"

using namespace matplot;

void calculate_thickness(std::vector<Point> points, int nx, int ny, float min_element_x, float min_element_y,
                         float step_length, float min_nx, float max_nx, float min_ny, float max_ny);

std::vector<Point> read_points_from_file(std::string path)
{

    std::string line;
    std::string value;

    std::vector<Point>
        points;
    std::ifstream file(path);

    if (file.is_open())
    {
        std::getline(file, line);
        while (std::getline(file, line))
        {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream linestream(line);
            Point p;
            linestream >> value;
            p.x = std::stod(value);
            linestream >> value;
            p.y = std::stod(value);
            linestream >> value;
            p.z = std::stod(value);
            // std::cout<< p.x << ","<< p.y<<"\n";

            points.emplace_back(p);
        }
    }
    else
    {
        std::cout << " File can't be opened \n";
    }
    file.close();

    return points;
}

std::vector<Width> sampling(std::vector<Point> points, double start, double x_end, double step_length)
{

    std::vector<Width> widthes;
    double start_x = start;
    double end = start + step_length;
    while (end < x_end)
    {
        int counter = 0;
        double min_y = 10000.;
        double max_y = -10000.;
        for (int i = 0; i < points.size(); i++)
        {
            double x = points[i].x;
            if (x >= start && x < end)
            {
                double y = points[i].y;
                if (y < min_y)
                    min_y = y;

                if (y > max_y)
                    max_y = y;

                counter++;
            }

            if (x > end)
                break;
        }

        // printf("%d, %f, %lf\n", counter, min_y, max_y);
        Width width;
        width.x = start + (step_length / 2.);
        width.left = max_y;
        width.right = min_y;
        width.length = max_y - min_y;

        widthes.emplace_back(width);

        start = end;
        end = start + step_length;
    }

    return widthes;
}

void save_to_csv(std::vector<Width> sampled_widthes, std::string path)
{

    std::ofstream file(path);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + path);
    }

    // Write the header row.
    file << "x,left,right,length" << std::endl;

    // Write the vector elements to the CSV file.
    for (auto element : sampled_widthes)
    {
        file << element.x << "," << element.left << "," << element.right << "," << element.length << std::endl;
    }

    file.close();
}

void plot_graph(std::vector<Point> points, std::vector<Width> widthes)
{

    std::vector<double> x;
    std::vector<double> y;

    std::vector<double> x_width;
    std::vector<double> left;
    std::vector<double> right;

    for (auto p : points)
    {
        x.push_back(p.x);
        y.push_back(p.y);
    }

    for (auto w : widthes)
    {
        x_width.push_back(w.x);
        left.push_back(w.left);
        right.push_back(w.right);
    }

    scatter(x, y, 15);
    hold(on);
    plot(x_width, left, "r-")->line_width(2);
    plot(x_width, right, "r-")->line_width(2);

    int counter = -1;
    double sum = 0;
    int count_sum = 0;
    float maxLength = 0.0;
    float minLength = 100000;
    for (auto w : widthes)
    {
        counter++;
        if (counter < 0.4 * widthes.size())
            continue;
        if (counter > 0.85 * widthes.size())
            continue;

        if (w.length > maxLength)
            maxLength = w.length;
        if (w.length < minLength)
            minLength = w.length;

        sum += w.length;
        count_sum++;
        std::vector<double> x{w.x, w.x};
        std::vector<double> y{w.left, w.right};
        int bold = 3;
        plot(x, y, "g-")->line_width(bold);
    }

    double avr = sum / float(count_sum);

    avr = avr + 0.75e-3;
    minLength = minLength + 0.75e-3;
    maxLength = maxLength + 0.75e-3;

    printf("min width:%lf mm\n", minLength * 1000.);
    printf("Average width:%lf mm\n", avr * 1000.);
    printf("Max width:%lf mm\n", maxLength * 1000.);

    show();
}

void plot_graph_3d(std::vector<Point> points, int nx, int ny, float min_element_x, float min_element_y, float step_length,
                   int min_nx, int max_nx,
                   int min_ny, int max_ny)
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> hash;

    for (auto p : points)
    {
        x.push_back(p.x);
        y.push_back(p.y);
        hash.push_back(p.hash);
    }

    scatter(x, y, 10, hash);
    hold(on);

    std::vector<float> x_grid;
    std::vector<float> y_grid;

    for (size_t i = min_nx; i <= max_nx; i++)
    {
        for (size_t j = min_ny; j <= max_ny; j++)
        {

            x_grid.push_back(min_element_x + i * step_length);
            y_grid.push_back(min_element_y + j * step_length);
        }
        plot(x_grid, y_grid, "r-")->line_width(1);
        x_grid.clear();
        y_grid.clear();
    }

    for (size_t j = min_ny; j <= max_ny; j++)
    {
        for (size_t i = min_nx; i <= max_nx; i++)
        {

            x_grid.push_back(min_element_x + i * step_length);
            y_grid.push_back(min_element_y + j * step_length);
        }
        plot(x_grid, y_grid, "r-")->line_width(1);
        x_grid.clear();
        y_grid.clear();
    }

    std::vector<float> thicknesses;
    for (size_t i = min_nx; i < max_nx; i++)
    {
        for (size_t j = min_ny; j < max_ny; j++)
        {

            float x_grid_min = min_element_x + i * step_length;
            float x_grid_max = x_grid_min + step_length;

            float y_grid_min = min_element_y + j * step_length;
            float y_grid_max = y_grid_min + step_length;

            float min_thickness = 1000;
            Point p_min;
            for (auto p : points)
            {
                if ((p.x > x_grid_min && p.x < x_grid_max) && (p.y > y_grid_min && p.y < y_grid_max))
                {
                    if (p.z < min_thickness)
                    {
                        min_thickness = p.z;
                        p_min = p;
                    }
                }
            }
            if (min_thickness != 1000)
                thicknesses.push_back(min_thickness);

            if (min_thickness * 1000. < -2.5)
            {
                std::vector<double> x;
                std::vector<double> y;
                x.push_back(p_min.x);
                y.push_back(p_min.y);
                scatter(x, y)->marker_face(true);
                x.clear();
                y.clear();
            }
        }
    }

    std::sort(thicknesses.begin(), thicknesses.end());

    std::cout << "******** Method (1) *********" << std::endl;

    float avr_thickness = std::reduce(thicknesses.begin(), thicknesses.end()) / thicknesses.size();
    std::cout << "avr_thickness: " << avr_thickness * -1000. << "\n";

    auto minElement = std::min_element(thicknesses.begin(), thicknesses.end());
    std::cout << "maxElement: " << *minElement * -1000. << "\n";

    auto maxElement = std::max_element(thicknesses.begin(), thicknesses.end());
    std::cout << "minElement: " << *maxElement * -1000. << "\n";

    calculate_thickness(points, nx, ny, min_element_x, min_element_y, step_length, min_nx, max_nx, min_ny, max_ny);

    show();
}
template <typename T>
void statistics(std::vector<T> &v)
{
    std::sort(v.begin(), v.end());
    /*     for (auto e : v)
        {
            std::cout << e << "\n";
        } */

    std::cout << "******** Method (2) *********" << std::endl;

    float avr_thickness = std::reduce(v.begin(), v.end()) / v.size();
    std::cout << "avr_thickness: " << avr_thickness * -1000. << "\n";

    auto minElement = std::min_element(v.begin(), v.end());
    std::cout << "maxElement: " << *minElement * -1000. << "\n";

    auto maxElement = std::max_element(v.begin(), v.end());
    std::cout << "minElement: " << *maxElement * -1000. << "\n";
}
void calculate_thickness(std::vector<Point> points, int nx, int ny, float min_element_x, float min_element_y, float step_length, float min_nx, float max_nx, float min_ny, float max_ny)
{
    std::vector<std::vector<float>> vec_thickness;
    for (size_t i = min_nx; i < max_nx; i++)
    {

        for (size_t j = min_ny; j < max_ny; j++)
        {

            float x_grid_min = min_element_x + i * step_length;
            float x_grid_max = x_grid_min + step_length;

            float y_grid_min = min_element_y + j * step_length;
            float y_grid_max = y_grid_min + step_length;

            std::vector<float> thicknesses;

            for (auto p : points)
            {
                if ((p.x > x_grid_min && p.x < x_grid_max) && (p.y > y_grid_min && p.y < y_grid_max))
                {
                    thicknesses.push_back(p.z);
                }
            }

            if (!thicknesses.empty())
                vec_thickness.push_back(thicknesses);
        }
    }

    std::vector<float> min_thickness;
    std::vector<int> count;
    for (auto vt : vec_thickness)
    {
        if (vt.size() < 10)
            continue;

        std::sort(vt.begin(), vt.end());
        auto minElement = (vt[0] + vt[1] + vt[2]) / 3.;
        min_thickness.push_back(minElement);
        count.push_back(vt.size());
    }

    statistics(min_thickness);
}

void compute_hashes(std::vector<Point> &points, long int nx, long int ny, long int nz, double dx,
                    double bbmin_x, double bbmin_y, double bbmin_z)
{

    for (size_t idx = 0; idx < points.size(); idx++)
    {
        float_t _px = points[idx].x;
        float_t _py = points[idx].y;
        float_t _pz = points[idx].z;

        int ix = (_px - bbmin_x) / dx;
        int iy = (_py - bbmin_y) / dx;
        int iz = (_pz - bbmin_z) / dx;

        // points[idx].hash = iz * ny * nx + iy * nx + ix;
        points[idx].hash = iy * nx + ix;
    }
}

/*  void unhash(int &i, int &j, int &k, int idx)
{

    k = idx / (geometry.ny * geometry.nx);
    j = (idx - k * geometry.ny * geometry.nx) / geometry.nx;
    i = idx % geometry.nx;
}
 */
int main()
{
    std::string fileName = "joined_particles_1500_rmft.csv";
    std::string path = "../Data/";
    path += fileName;

    std::cout << fileName << std::endl;
    std::vector<Point> points = read_points_from_file(path);

    /* std::cout << points.size() << "\n";
    for (auto p : points)
        std::cout << p.x << "," << p.y << "\n"; */

    auto min_element_iterator_x = std::min_element(points.begin(), points.end(), compare_x);
    auto max_element_iterator_x = std::max_element(points.begin(), points.end(), compare_x);

    auto min_element_iterator_y = std::min_element(points.begin(), points.end(), compare_y);
    auto max_element_iterator_y = std::max_element(points.begin(), points.end(), compare_y);

    auto min_element_iterator_z = std::min_element(points.begin(), points.end(), compare_z);
    auto max_element_iterator_z = std::max_element(points.begin(), points.end(), compare_z);

    double min_element_x = (*min_element_iterator_x).x;
    double max_element_x = (*max_element_iterator_x).x;

    double min_element_y = (*min_element_iterator_y).y;
    double max_element_y = (*max_element_iterator_y).y;

    double min_element_z = (*min_element_iterator_z).z;
    double max_element_z = (*max_element_iterator_z).z;

    // printf("%f, %f, %f, %f, %f, %f \n", min_element_x, max_element_x, min_element_y, max_element_y, min_element_z, max_element_z);

    double x_distance = max_element_x - min_element_x;

    // double sample = 20.;
    // double step_length = x_distance / sample;
    double step_length = 0.002;

    printf("step_length:  %lf \n", step_length);

    std::sort(points.begin(), points.end(), compare_x);
    std::vector<Width> sampled_widthes = sampling(points, min_element_x, max_element_x, step_length);
    save_to_csv(sampled_widthes, "../out/sampled_width.csv");

    // plot_graph(points, sampled_widthes);

    int nx = ceil(std::abs(max_element_x - min_element_x) / step_length);
    int ny = ceil(std::abs(max_element_y - min_element_y) / step_length);
    int nz = 1;

    // printf("Number of cells: %d, %d, %d, %lf \n", nx, ny, nz, step_length);
    printf("Number of cells: %d, %d, %d, %lf \n", nx, ny, nz, step_length);

    compute_hashes(points, nx, ny, nz, step_length,
                   min_element_x, min_element_y, min_element_z);

    std::sort(points.begin(), points.end(), compare_hash);

    // start
    int min_nx = 0.0 * nx;
    int max_nx = nx;

    int min_ny = 0 * ny;
    int max_ny = 1. * ny;

    plot_graph_3d(points, nx, ny, min_element_x, min_element_y, step_length, min_nx, max_nx, min_ny, max_ny);
    // calculate_thickness(points, nx, ny, min_element_x, min_element_y, step_length);
}
