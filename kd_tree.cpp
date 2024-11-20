#include "kd_tree.h"
#include <algorithm>
#include <iostream>
#include <cassert>

template <typename T>
void kd_tree<T>::build_tree(const std::vector<std::vector<T>>& data)
{
    std::vector<kd_point<T>> copy;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
        copy.emplace_back(i, data[i]);
    }
    root_ = std::make_unique<kd_node<T>>();
    root_->build_node(copy, copy.begin(), copy.end());
}

template <typename T>
kd_point<T> kd_tree<T>::closest_neighbour(const std::vector<T>& point_query)
{
    if (root_ == nullptr) return kd_point<T>();
    kd_point<T> closest;
    closest.idx_ = limit<int>::min();
    closest.point_coords_.resize(point_query.size());
    T max = limit<T>::max();
    root_->closest_neighbour(point_query, closest, max);
    return closest;
}

template <typename T>
void kd_node<T>::build_node(std::vector<kd_point<T>>& data, kd_iterator<T> start, kd_iterator<T> end, int separation_axis)
{
    const int64_t size = end - start;
    separation_axis_ = (separation_axis + 1) % start->point_coords_.size();
    if (size < 2)
    {
        point_ = *start;
        left_ = nullptr;
        right_ = nullptr;
    }
    else
    {
        kd_iterator<T> median = start + size / 2;
        std::nth_element(start, median, end,
            [this](const kd_point<T>& point1, const kd_point<T>& point2) -> bool
            {
                return point1.point_coords_[separation_axis_] < point2.point_coords_[separation_axis_];
            });
        point_ = *median;
        if (start != median)
        {
            left_ = std::make_unique<kd_node>();
            left_->build_node(data, start, median, separation_axis_);
        }
        if (median != end)
        {
            right_ = std::make_unique<kd_node>();
            right_->build_node(data, median, end, separation_axis_);
        }
    }
}

template <typename T>
void kd_node<T>::closest_neighbour(const std::vector<T>& point_query, kd_point<T>& closest, double& closest_sqr_euclidean_dist)
{
    assert(point_.point_coords_.size() == point_query.size());
    T dist = helper_func::squared_euclidean_distance<T>(point_.point_coords_, point_query);
    if (dist < closest_sqr_euclidean_dist )
    {
        closest = point_;
        closest_sqr_euclidean_dist = dist;
    }
    if (left_ != nullptr && point_query[separation_axis_] <= point_.point_coords_[separation_axis_])
        left_->closest_neighbour(point_query, closest, closest_sqr_euclidean_dist);
    else if(right_ != nullptr)
        right_->closest_neighbour(point_query, closest, closest_sqr_euclidean_dist);
    // check the other node when backtracking
    const double other_side_dist = point_query[separation_axis_] - point_.point_coords_[separation_axis_];
    if (other_side_dist * other_side_dist <= closest_sqr_euclidean_dist)
    {
        if (right_ != nullptr && point_query[separation_axis_] <= point_.point_coords_[separation_axis_])
            right_->closest_neighbour(point_query, closest, closest_sqr_euclidean_dist);
        else if (left_ != nullptr)
            left_->closest_neighbour(point_query, closest, closest_sqr_euclidean_dist);
    }
}

template<typename T>
double helper_func::squared_euclidean_distance(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());
    double distance = 0.0;
    for (unsigned int i = 0; i < a.size(); ++i)
    {
        distance += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return distance;
}

template<typename T>
void helper_func::print_example(const std::vector<std::vector<T>>& list, const std::vector<T>& query, const kd_point<T>& point)
{
    std::cout << "Tree constructed from points: \n[";
    for(unsigned int i = 0; i < list.size(); ++i)
    {
        std::cout << "(";
        for(unsigned int j = 0; j < list[i].size(); ++j)
        {
            j != list[i].size() -1 ? (std::cout << list[i][j] << ",") : (std::cout << list[i][j]);
        }
        i != list.size() - 1 ? std::cout << ")," : std::cout << ")";
    }
    std::cout << "]\n";
    std::cout << "Query: \n(";
    for(unsigned int i = 0; i < query.size(); ++i)
    {
        i != query.size() -1 ? (std::cout << query[i] << ",") : (std::cout << query[i]);
    }
    std::cout << ")\n";
    std::cout << "Closest neighbour to searched value is:\n(";
    for(unsigned int i = 0; i < point.point_coords_.size(); ++i)
    {
        i != point.point_coords_.size() -1 ? (std::cout << point.point_coords_[i] << ",") : (std::cout << point.point_coords_[i]);
    }
    std::cout << ")\n";
}

int main(int argc, char* argv[])
{
    kd_tree<double> tree;
    std::vector<std::vector<double>> test {{9,1}, {6,12}, {3,6}, {13,15}, {10,19}, {17,15}, {15,12}};
    tree.build_tree(test);
    std::vector<double> query {2, 5};
    auto res = tree.closest_neighbour(query);
    helper_func::print_example(test, query, res);
    test = {{9,1,5}, {6,12,1}, {3,6,4}, {13,15,2}, {10,19,0}, {17,15,6}, {15,12,-12}};
    tree.build_tree(test);
    query = {2, 5, -7};
    res = tree.closest_neighbour(query);
    helper_func::print_example(test, query, res);
    return 0;
}
