#pragma once
#include <vector>
#include <memory>

template<typename T>
using limit = std::numeric_limits<T>;

template<typename T>
class kd_point;

template<typename T>
using kd_iterator = typename std::vector<kd_point<T>>::iterator;

template<typename T>
class kd_point
{
public:
    kd_point() : idx_(limit<int>::min()) {}
    kd_point(int idx, const std::vector<T>& point_coords) : idx_(idx), point_coords_(point_coords) {}

    int idx_;
    std::vector<T> point_coords_;
};

template<typename T>
class kd_node
{
public:
    kd_node() : left_(nullptr), right_(nullptr) {}
    void build_node(std::vector<kd_point<T>>& data,
       kd_iterator<T> start,
       kd_iterator<T> end,
       int separation_axis = -1
   );
    void closest_neighbour(const std::vector<T>& point_query,
        kd_point<T>& closest, double& closest_sqr_euclidean_dist);
private:
    std::unique_ptr<kd_node> left_;
    std::unique_ptr<kd_node> right_;
    kd_point<T> point_;
    int separation_axis_{-1};
};

template<typename T>
class kd_tree
{
public:
    kd_tree() : root_(nullptr) {}
    void build_tree(const std::vector<std::vector<T> >& data);
    kd_point<T> closest_neighbour(const std::vector<T>& point_query);
private:
    std::unique_ptr<kd_node<T>> root_;
};

namespace helper_func
{
    template<typename T>
    double squared_euclidean_distance(const std::vector<T>& a, const std::vector<T>& b);

    template<typename T>
    void print_example(const std::vector<std::vector<T>>& list, const std::vector<T>& query, const kd_point<T>& point);
};

