#pragma once

#include <string_view>
#include <fstream>

#include "Geometry.cpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

class OFFReader {
public:
    static Geometry read(const std::string_view filename) {

        std::ifstream off_file(filename.data());

        std::string off_file_start_token = "OFF";
        unsigned int num_vertices = 0;
        unsigned int num_faces = 0;
        unsigned int num_edges = 0;

        //Skip the line with the token that signals the file start ("OFF")
        std::getline(off_file, off_file_start_token);

        off_file >> num_vertices >> num_faces >> num_edges;

        std::vector<Kernel::Point_2> vertices;
        std::vector<std::vector<unsigned int>> faces;
        std::vector<std::pair<unsigned int, unsigned int>> edges;

        vertices.reserve(num_vertices);
        faces.reserve(num_faces);

        for(int vertex_index = 0; vertex_index != num_vertices; ++vertex_index){
            double x = 0.0;
            double y = 0.0;

            off_file >> x >> y;

            vertices.emplace_back(x, y);
        }

        for(int face_index = 0; face_index != num_faces; ++face_index){
            unsigned int num_vertices_in_face = 0;

            off_file >> num_vertices_in_face;

            auto& last_vector_added = faces.emplace_back();
            last_vector_added.reserve(num_vertices_in_face);

            for(int i = 0; i != num_vertices_in_face; ++i){
                unsigned int face_vertex_index = 0;
                off_file >> face_vertex_index;
                last_vector_added.push_back(face_vertex_index);
            }
        }

        for(const auto& face_vertices : faces){
            for(int i = 0; i < face_vertices.size(); ++i){
                edges.emplace_back(std::pair(face_vertices[i], face_vertices[(i + 1) % face_vertices.size()]));
            }
        }

        return Geometry(std::move(vertices), std::move(faces), std::move(edges));
    }
    static void write_dcel(const std::string_view filename, const Arrangement_2& dcel) {

        std::ofstream output_off_file(filename.data());

        std::string off_file_start_token = "OFF\n";
        unsigned int num_vertices = dcel.number_of_vertices();
        unsigned int num_faces = dcel.number_of_faces();
        unsigned int num_edges = dcel.number_of_edges();

        //Skip the line with the token that signals the file start ("OFF")
        output_off_file << off_file_start_token;

        output_off_file << num_vertices << " " << num_faces << " " << num_edges << "\n";

    }

};        