#include <iostream>
#include <fstream>
#include <queue>

#define DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "CGAL/Point_2.h"
#include "CGAL/Arrangement_on_surface_2.h"
#include "CGAL/Arr_segment_traits_2.h"
#include "CGAL/Arrangement_2.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
using Vertex_handle = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Vertex_handle;
using Vertex = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Vertex;
using Half_edge = CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>::Halfedge;


typedef int                                           Number_type;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef CGAL::Point_2<Kernel>                         Point_2;
typedef CGAL::Line_2<Kernel>                          Line_2;
typedef CGAL::Direction_2<Kernel>                     Direction_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Face_handle                    Face_handle;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;
typedef Arrangement_2::Vertex_const_handle            Vertex_const_handle;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef Arrangement_2::Halfedge_const_handle          Halfedge_const_handle;

enum class event_type {
    start,
    end,
    upper_regular,
    lower_regular,
    split,
    merge
};

void print_edge(Halfedge_handle edge){
    std::cerr << edge->source()->point().x() << " " << edge->source()->point().y() << " -> ";
    std::cerr << edge->target()->point().x() << " " << edge->target()->point().y() << "  \n";
}
void print_vertex(Vertex_handle v){
    std::cerr << v->point().x() << " " << v->point().y() << "  \n";
}
void print_point(Point_2 v){
    std::cerr << v.x() << " " << v.y() << "  \n";
}

struct vertex_event {
    event_type type;
    Vertex_handle vertex_handle;
    Halfedge_handle next_edge;
    Halfedge_handle prev_edge;
};

class EdgeStructure {
private:
    class EdgeOrder {

        Line_2& sl;

    public:

        EdgeOrder(Line_2& _sl) : sl(_sl) {}

        bool operator() (Halfedge_handle i, Halfedge_handle j ) const {

            const Point_2 *p = nullptr;
            const Point_2 *q = nullptr;

            CGAL::Object u = CGAL::intersection( CGAL::Line_2<Kernel>(i->source()->point(), i->target()->point()), sl );
            CGAL::Object v = CGAL::intersection( CGAL::Line_2<Kernel>(j->source()->point(), j->target()->point()), sl );

            if ( CGAL::object_cast<Point_2>(&u) ) {
                p = CGAL::object_cast<Point_2>( &u );
            } else if ( CGAL::object_cast<Line_2>(&u) ) {
                p = (const Point_2*) i->source().ptr();
            }
            if ( CGAL::object_cast<Point_2>(&v) ) {
                q = CGAL::object_cast<Point_2>( &v );
            } else if ( CGAL::object_cast<Line_2>(&v) ) {
                q = (const Point_2*) j->source().ptr();
            }
            if ( (p != nullptr) && (q != nullptr) ) {

                return ( (*p).y() < (*q).y() );

            } else {

                return false;
            }
        }

    };

    Line_2 sl;
    std::multiset<Halfedge_handle, EdgeOrder> edges;
public:

    EdgeStructure() : edges((EdgeOrder(sl))) {}

    bool empty() const {
        return edges.empty();
    }

    void insert(const Halfedge_handle& e){
        sl = Line_2(e->source()->point(), Direction_2(0, 1));
        edges.insert(e);
    }
    void erase(const Halfedge_handle& e){
        edges.erase(e);
    }

    //Troviamo l'edge direttamente sotto il vertice
    Halfedge_handle find_edge_underneath(const Halfedge_handle& e){
        if(edges.size() == 1) return *edges.begin();
        sl = Line_2(e->source()->point(), Direction_2(0, 1));
        auto it = edges.upper_bound(e);
        return *std::prev(it);
    }

    friend std::ostream& operator<<(std::ostream& out, const EdgeStructure& edge_structure){
        out << "-----------------\n";
        for(const auto e : edge_structure.edges){
            out << e->source()->point().x() << " " << e->source()->point().y();
            out << " -> ";
            out << e->target()->point().x() << " " << e->target()->point().y();
            out << "\n";
        }
        return out;
    }

};

void print_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    std::cout << "Vertices\n";
    for(Vertex_handle& v : dcel.vertex_handles()){
        std::cout << v->point().x() << " " << v->point().y() << "\n";
    }

    int face_index = 0;
    for(const auto& face : dcel.face_handles()){
        if(face->is_unbounded()) {
            face_index++;
            continue;
        }

        std::cout << "Face: " << face_index++ << "\n";
        auto edge_circulator = face->outer_ccb();
        do {
            auto he = edge_circulator;
            std::cout << he->source()->point().x() << " " << he->source()->point().y() << " -> " << he->target()->point().x() << " " << he->target()->point().y() << "\n";

        } while (++edge_circulator != face->outer_ccb());
    }

    std::cout << "Edges\n";
    for(const auto& he : dcel.edge_handles()){
        std::cout << he->source()->point().x() << " " << he->source()->point().y() << " -> " << he->target()->point().x() << " " << he->target()->point().y() << "\n";
    }

}

void insert(const Halfedge_handle& e,
            EdgeStructure& sweep_line_edges,
            std::map<Halfedge_const_handle, Vertex_handle>& helpers,
            const std::map<Vertex_const_handle, event_type>& types){
    sweep_line_edges.insert(e);
    helpers[e] = e->source();
}
void remove(const Halfedge_handle& e,
            const Vertex_handle& v,
            EdgeStructure& sweep_line_edges,
            std::map<Halfedge_const_handle, Vertex_handle>& helpers,
            const std::map<Vertex_const_handle, event_type>& types,
            CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){
    const auto& u = helpers[e];

    if(types.contains(u) && types.at(u) == event_type::merge){
        //Aggiungiamo la diagonale da u a v
        const Segment_2 segment(u->point(), v->point());
#ifdef DEBUG
        std::cerr << "Insert diagonal ";
        std::cerr << u->point().x() << " " << u->point().y() << " -> ";
        std::cerr << v->point().x() << " " << v->point().y() << "";
        std::cerr << "\n";
#endif
        dcel.insert_at_vertices(segment, u, v);
    }

    sweep_line_edges.erase(e);
}
void process(const Vertex_handle& v,
             const Halfedge_handle& e,
             EdgeStructure& sweep_line_edges,
             std::map<Halfedge_const_handle, Vertex_handle>& helpers,
             const std::map<Vertex_const_handle, event_type>& types,
             CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    const auto b = sweep_line_edges.find_edge_underneath(e);
#ifdef DEBUG
    std::cerr << "Found edge ";
    print_edge(b);
#endif
    const auto& u = helpers.at(b);
    if(types.at(u) == event_type::merge || types.at(v) == event_type::split){
        //Aggiunta della diagonale uv
        const Segment_2 segment(u->point(), v->point());
#ifdef DEBUG
        std::cerr << "Insert diagonal ";
        std::cerr << u->point().x() << " " << u->point().y() << " -> ";
        std::cerr << v->point().x() << " " << v->point().y() << "";
        std::cerr << "\n";
#endif
        dcel.insert_at_vertices(segment, u, v);
    }

    helpers[b] = v;
}

event_type identify_event(const Point_2& prev_point, const Point_2& point, const Point_2& next_point){

    if (CGAL::lexicographically_xy_smaller(prev_point, point) &&
        CGAL::lexicographically_xy_smaller(point, next_point)) {
        //lower regular
        return event_type::upper_regular;
    }
    if (CGAL::lexicographically_xy_larger(prev_point, point) &&
        CGAL::lexicographically_xy_larger(point, next_point)) {
        //upper regular
        return event_type::lower_regular;
    }
    if (CGAL::lexicographically_xy_larger(prev_point, point) &&
        CGAL::lexicographically_xy_larger(next_point, point)) {
        //Based on turn direction identify start or split event
        if (CGAL::left_turn(prev_point, point, next_point)) {
            return event_type::split;
        } else {
            return event_type::start;
        }
    }
    if (CGAL::lexicographically_xy_smaller(prev_point, point) &&
        CGAL::lexicographically_xy_smaller(next_point, point)) {
        //Based on turn direction identify end or merge event
        if (CGAL::left_turn(prev_point, point, next_point)) {
            return event_type::merge;
        } else {
            return event_type::end;
        }
    }

}

void monotone_subdivision(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel, Face_handle f){

    //For every vertex register the event types
    std::map<Vertex_const_handle, event_type> types;

    std::vector<vertex_event> vertex_events;
    vertex_events.reserve(dcel.number_of_vertices());

    auto current = f->outer_ccb();
    std::cout << "Next face\n";
    do {
        auto edge_it = current;

        Halfedge_handle prev_handle = (++edge_it)->twin();
        Halfedge_handle curr_handle = current->twin();

        auto prev = *prev_handle;
        auto edge = *curr_handle;

        auto vertex = edge.source();
        const auto next_vertex = edge.target();
        const auto prev_vertex = prev.source();
        const auto point = vertex->point();
        const auto next_point = next_vertex->point();
        const auto prev_point = prev_vertex->point();

        const auto event_type = identify_event(prev_point, point, next_point);

        vertex_events.emplace_back(event_type, vertex, prev_handle->twin(), curr_handle->twin());
#ifdef DEBUG
        std::cerr << "Outer edge \n";
        print_vertex(vertex);
        print_edge(prev_handle->twin());
        print_edge(curr_handle->twin());
#endif
        types[vertex] = event_type;
    } while (++current != f->outer_ccb());

    for (auto hi = f->holes_begin(); hi != f->holes_end(); ++hi) {
        auto current_inner = *hi;
        do {
            auto inverse_handle = current_inner;
            auto prev_handle = inverse_handle;

            auto edge = *(inverse_handle);

            auto vertex = edge.source();
            const auto next_vertex = edge.target();
            const auto prev_vertex = (--prev_handle)->source();

            const auto point = vertex->point();
            const auto next_point = next_vertex->point();
            const auto prev_point = prev_vertex->point();

            const auto event_type = identify_event(next_point, point, prev_point);

            vertex_events.emplace_back(
                    event_type,
                    vertex,
                    inverse_handle,
                    inverse_handle->prev()
            );
            #ifdef DEBUG
            std::cerr << "Inner edges \n";
            print_vertex(vertex);
            print_edge(inverse_handle);
            print_edge(inverse_handle->prev());
            #endif
            types[vertex] = event_type;
        } while (++current_inner != *hi);
    }

    //Ordiniamo gli eventi da sinistra a destra
    std::sort(vertex_events.begin(), vertex_events.end(), [](const auto& a, const auto& b){
        return CGAL::lexicographically_xy_smaller(a.vertex_handle->point(), b.vertex_handle->point());
    });

    std::map<Halfedge_const_handle, Vertex_handle> helpers;
    EdgeStructure sweep_line_edges;

    //Processiamo ciascun evento
    for(const auto& event : vertex_events) {

        const auto event_type = event.type;

#ifdef DEBUG
        print_vertex(event.vertex_handle);
#endif

        if (event_type == event_type::start) {
#ifdef DEBUG
            std::cerr << "event_type::start\n";
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
#endif
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::end) {
#ifdef DEBUG
            std::cerr << "event_type::end\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
#endif
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
        }
        if (event_type == event_type::lower_regular) {
#ifdef DEBUG
            std::cerr << "event_type::lower_regular\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
#endif
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::upper_regular) {
#ifdef DEBUG
            std::cerr << "event_type::upper_regular\n";
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
#endif
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }
        if (event_type == event_type::split) {
#ifdef DEBUG
            std::cerr << "event_type::split\n";
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
#endif
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::merge) {
#ifdef DEBUG
            std::cerr << "event_type::merge\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
#endif
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
#ifdef DEBUG
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
#endif
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }
    }

}

CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> example_1(){
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    auto v1 = dcel.insert_in_face_interior(Point_2(1.0, 2.0), dcel.unbounded_face());
    auto v2 = dcel.insert_in_face_interior(Point_2(3.0, 1.0), dcel.unbounded_face());
    auto v3 = dcel.insert_in_face_interior(Point_2(2.0, 2.0), dcel.unbounded_face());
    auto v4 = dcel.insert_in_face_interior(Point_2(3.0, 3.0), dcel.unbounded_face());

    dcel.insert_at_vertices(Segment_2(v1->point(), v2->point()), v1, v2);
    dcel.insert_at_vertices(Segment_2(v2->point(), v3->point()), v2, v3);
    dcel.insert_at_vertices(Segment_2(v3->point(), v4->point()), v3, v4);
    dcel.insert_at_vertices(Segment_2(v4->point(), v1->point()), v4, v1);

    return dcel;
}
CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> example_2(){
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    auto v1 = dcel.insert_in_face_interior(Point_2(1.0, 6.0), dcel.unbounded_face());
    auto v2 = dcel.insert_in_face_interior(Point_2(3.0, 1.0), dcel.unbounded_face());
    auto v3 = dcel.insert_in_face_interior(Point_2(8.0, 1.0), dcel.unbounded_face());
    auto v4 = dcel.insert_in_face_interior(Point_2(10.0, 4.0), dcel.unbounded_face());
    auto v5 = dcel.insert_in_face_interior(Point_2(12.0, 1.0), dcel.unbounded_face());
    auto v6 = dcel.insert_in_face_interior(Point_2(14.0, 1.0), dcel.unbounded_face());
    auto v7 = dcel.insert_in_face_interior(Point_2(10.0, 7.0), dcel.unbounded_face());
    auto v8 = dcel.insert_in_face_interior(Point_2(12.0, 7.0), dcel.unbounded_face());
    auto v9 = dcel.insert_in_face_interior(Point_2(13.0, 5.0), dcel.unbounded_face());
    auto v10 = dcel.insert_in_face_interior(Point_2(13.0, 8.0), dcel.unbounded_face());
    auto v11 = dcel.insert_in_face_interior(Point_2(4.0, 8.0), dcel.unbounded_face());

    dcel.insert_at_vertices(Segment_2(v1->point(), v2->point()), v1, v2);
    dcel.insert_at_vertices(Segment_2(v2->point(), v3->point()), v2, v3);
    dcel.insert_at_vertices(Segment_2(v3->point(), v4->point()), v3, v4);
    dcel.insert_at_vertices(Segment_2(v4->point(), v5->point()), v4, v5);
    dcel.insert_at_vertices(Segment_2(v5->point(), v6->point()), v5, v6);
    dcel.insert_at_vertices(Segment_2(v6->point(), v7->point()), v6, v7);
    dcel.insert_at_vertices(Segment_2(v7->point(), v8->point()), v7, v8);
    dcel.insert_at_vertices(Segment_2(v8->point(), v9->point()), v8, v9);
    dcel.insert_at_vertices(Segment_2(v9->point(), v10->point()), v9, v10);
    dcel.insert_at_vertices(Segment_2(v10->point(),v11->point()), v10, v11);
    dcel.insert_at_vertices(Segment_2(v11->point(),v1->point()), v11, v1);

    return dcel;
}
CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> example_3(){
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    std::vector<Point_2> boundary_points = {
            Point_2(0, 5),
            Point_2(5, 0),
            Point_2(10, 5),
            Point_2(5, 10),
    };

    std::vector<Point_2> internal_points = {
            Point_2(2.5, 5),
            Point_2(5, 2.5),
            Point_2(7.5, 5),
            Point_2(5, 7.5),
    };

    Vertex_handle vh0 = dcel.insert_in_face_interior(boundary_points[0], dcel.unbounded_face());
    Vertex_handle vh1 = dcel.insert_in_face_interior(boundary_points[1], dcel.unbounded_face());
    Vertex_handle vh2 = dcel.insert_in_face_interior(boundary_points[2], dcel.unbounded_face());
    Vertex_handle vh3 = dcel.insert_in_face_interior(boundary_points[3], dcel.unbounded_face());

    Vertex_handle ivh0 = dcel.insert_in_face_interior(internal_points[0], dcel.faces_begin());
    Vertex_handle ivh1 = dcel.insert_in_face_interior(internal_points[1], dcel.faces_begin());
    Vertex_handle ivh2 = dcel.insert_in_face_interior(internal_points[2], dcel.faces_begin());
    Vertex_handle ivh3 = dcel.insert_in_face_interior(internal_points[3], dcel.faces_begin());

    auto h1 = dcel.insert_at_vertices(Segment_2(boundary_points[0], boundary_points[1]), vh0, vh1);
    auto h2 = dcel.insert_at_vertices(Segment_2(boundary_points[1], boundary_points[2]), h1, vh2);
    auto h3 = dcel.insert_at_vertices(Segment_2(boundary_points[2], boundary_points[3]), h2, vh3);
    auto h4 = dcel.insert_at_vertices(Segment_2(boundary_points[3], boundary_points[0]), h3, vh0);

    auto ih1 = dcel.insert_at_vertices(Segment_2(internal_points[0], internal_points[3]), ivh0, ivh3);
    dcel.insert_at_vertices(Segment_2(internal_points[3], internal_points[2]), ivh3, ivh2);
    dcel.insert_at_vertices(Segment_2(internal_points[2], internal_points[1]), ivh2, ivh1);
    dcel.insert_at_vertices(Segment_2(internal_points[0], internal_points[1]), ivh0, ivh1);

    return dcel;
}

CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> example_4(){
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    std::array<Point_2, 38> points = {
        Point_2(0.5, 4.5),
        Point_2(1.5, 2),
        Point_2(4.5, 0.5),
        Point_2(7.25, 1.75),
        Point_2(6.25, 2.10),
        Point_2(7.5, 3),
        Point_2(8.25, 0.5),
        Point_2(11.5, 1),
        Point_2(12.6, 4),
        Point_2(11.5, 5.5),
        Point_2(7.5, 4.6),
        Point_2(9.5, 6),
        Point_2(5.75, 5.5),
        Point_2(3.5, 6),

        Point_2(1.75, 3.5),
        Point_2(3, 1.525),

        Point_2(2.6, 3.8),
        Point_2(3.10, 2.75),
        Point_2(4.5, 1.5),
        Point_2(6, 2),
        Point_2(5.5, 2.6),
        Point_2(6.5, 3.75),
        Point_2(8.25, 3.25),
        Point_2(9.75, 2.15),
        Point_2(8.8, 3.25),
        Point_2(9.75, 4.6),
        Point_2(8.25, 4),
        Point_2(5.75, 4.75),

        Point_2(10.40, 1.8),
        Point_2(11.40, 2.8),
        Point_2(10.75, 3.25),
        Point_2(11.50, 4),
        Point_2(10.60, 4.5),
        Point_2(10.10, 3.25),

        Point_2(9.5, 3.25),

        Point_2(4.20, 4),
        Point_2(4.5, 2.25),
        Point_2(6, 3.8)

    };

    std::array<Vertex_handle , points.size()> vertices;
    for(int i = 0; i < 16; i++) vertices[i] = dcel.insert_in_face_interior(points[i], dcel.unbounded_face());

    //Faccia 1
    dcel.insert_at_vertices(Segment_2(vertices[0]->point(), vertices[1]->point()), vertices[0], vertices[1]);
    dcel.insert_at_vertices(Segment_2(vertices[1]->point(), vertices[2]->point()), vertices[1], vertices[2]);
    dcel.insert_at_vertices(Segment_2(vertices[2]->point(), vertices[15]->point()), vertices[2], vertices[15]);
    dcel.insert_at_vertices(Segment_2(vertices[15]->point(), vertices[14]->point()), vertices[15], vertices[14]);
    dcel.insert_at_vertices(Segment_2(vertices[14]->point(), vertices[13]->point()), vertices[14], vertices[13]);
    dcel.insert_at_vertices(Segment_2(vertices[13]->point(), vertices[0]->point()), vertices[13], vertices[0]);

    //Faccia 2
    for(int i = 2; i < 13; i++) dcel.insert_at_vertices(Segment_2(vertices[i]->point(), vertices[i+1]->point()), vertices[i], vertices[i+1]);

    auto right_face = ++(++dcel.faces_begin());
    for(int i = 16; i <= 34; i++) vertices[i] = dcel.insert_in_face_interior(points[i], right_face);

    //Faccia 3
    for(int i = 16; i < 27; i++) dcel.insert_at_vertices(Segment_2(vertices[i]->point(), vertices[i+1]->point()), vertices[i], vertices[i+1]);
    dcel.insert_at_vertices(Segment_2(vertices[27]->point(), vertices[16]->point()), vertices[27], vertices[16]);

    //Faccia 4 e 5
    for(int i = 28; i < 33; i++) dcel.insert_at_vertices(Segment_2(vertices[i]->point(), vertices[i+1]->point()), vertices[i], vertices[i+1]);
    dcel.insert_at_vertices(Segment_2(vertices[33]->point(), vertices[28]->point()), vertices[33], vertices[28]);

    dcel.insert_at_vertices(Segment_2(vertices[32]->point(), vertices[34]->point()), vertices[32], vertices[34]);
    dcel.insert_at_vertices(Segment_2(vertices[34]->point(), vertices[28]->point()), vertices[34], vertices[28]);

    auto third_face = (++(++(++dcel.faces_begin())));
    for(int i = 35; i < points.size(); i++) vertices[i] = dcel.insert_in_face_interior(points[i], third_face);
    dcel.insert_at_vertices(Segment_2(vertices[35]->point(), vertices[36]->point()), vertices[35], vertices[36]);
    dcel.insert_at_vertices(Segment_2(vertices[36]->point(), vertices[37]->point()), vertices[36], vertices[37]);
    dcel.insert_at_vertices(Segment_2(vertices[37]->point(), vertices[35]->point()), vertices[37], vertices[35]);

    return dcel;
}

int main() {

    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel_1 = example_1();
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel_2 = example_2();
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel_3 = example_3();
    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel_4 = example_4();

    print_dcel(dcel_4);
    //monotone_subdivision(dcel_1, ++dcel_1.faces_begin()); //example 1
    //monotone_subdivision(dcel_2, ++dcel_2.faces_begin()); //example 2
    //monotone_subdivision(dcel_3, ++dcel_3.faces_begin()); //example 3
    auto f = dcel_4.faces_begin();
    while(f != dcel_4.faces_end()) {
        if (f->is_unbounded()) {
            ++f;
            continue;
        }
        monotone_subdivision(dcel_4, f); //example 4
        ++f;
    }
    print_dcel(dcel_4);

    return 0;
}