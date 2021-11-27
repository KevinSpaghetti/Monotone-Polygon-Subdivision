//
// Created by Kevin on 17/11/2021.
//

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "CGAL/Point_2.h"
#include "CGAL/Arr_segment_traits_2.h"
#include "CGAL/Arrangement_2.h"
#include "CGAL/Surface_mesh.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
using Vertex_handle = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Vertex_handle;
using Vertex = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Vertex;
using Half_edge = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Halfedge;

typedef int                                           Number_type;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef CGAL::Point_2<Kernel>                         Point_2;
typedef CGAL::Line_2<Kernel>                          Line_2;
typedef CGAL::Direction_2<Kernel>                          Direction_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Vertex_const_handle            Vertex_const_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef Arrangement_2::Halfedge_const_handle          Halfedge_const_handle;

struct Geometry{
    std::vector<Point_2> vertices;
    std::vector<std::vector<unsigned int>> faces;
    std::vector<std::pair<unsigned int, unsigned int>> edges;

    Geometry(
            std::vector<Point_2>&& _vertices,
            std::vector<std::vector<unsigned int>>&& _faces,
            std::vector<std::pair<unsigned int, unsigned int>>&& _edges) :
            vertices(std::forward<std::vector<Point_2>>(_vertices)),
            faces(std::forward<std::vector<std::vector<unsigned int>>>(_faces)),
            edges(std::forward<std::vector<std::pair<unsigned int, unsigned int>>>(_edges)){}

    Geometry(
            const std::vector<Point_2>& _vertices,
            const std::vector<std::vector<unsigned int>>& _faces,
            const std::vector<std::pair<unsigned int, unsigned int>>& _edges) :
            vertices(_vertices),
            faces(_faces),
            edges(_edges){}

};