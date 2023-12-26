struct Point

{
    double x;
    double y;
    double z;
    int hash;
};

struct Width
{
    double x;
    double left;
    double right;
    double length;
};

bool compare_x(Point a, Point b)
{
    return a.x < b.x;
}

bool compare_y(Point a, Point b)
{
    return a.y < b.y;
}

bool compare_z(Point a, Point b)
{
    return a.z < b.z;
}

bool compare_hash(Point a, Point b)
{
    return a.hash< b.hash;
}