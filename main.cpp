#include <iostream>
#include <fstream>
#include <queue>

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

struct vertex_event {
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

    void insert(const Halfedge_handle& e){
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

//Stampa la DCEL in un formato leggibile
void print_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){


    std::cout << "Vertici\n";
    for(Vertex_handle& v : dcel.vertex_handles()){
        std::cout << v->point().x() << " " << v->point().y() << "\n";
    }

    int face_index = 0;
    for(const auto& face : dcel.face_handles()){
        if(face->is_unbounded()) {
            face_index++;
            continue;
        }

        std::cout << "Faccia: " << face_index++ << "\n";
        auto edge_circulator = face->outer_ccb();
        do {
            auto he = edge_circulator;
            std::cout << he->source()->point().x() << " " << he->source()->point().y() << " -> " << he->target()->point().x() << " " << he->target()->point().y() << "\n";

        } while (++edge_circulator != face->outer_ccb());
    }

    std::cout << "Half-edges\n";
    for(const auto& he : dcel.edge_handles()){
        std::cout << he->source()->point().x() << " " << he->source()->point().y() << " -> " << he->target()->point().x() << " " << he->target()->point().y() << "\n";
    }

}

//Leggere e stampare la DCEL utilizzando lo standard input (in formato OFF)
CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> read_dcel(){

    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    std::string off_file_start_token = "OFF";
    unsigned int num_vertices = 0;
    unsigned int num_faces = 0;
    unsigned int num_edges = 0;

    //Saltiamo il token che segnala l'inizio di un file OFF ("OFF")
    std::getline(std::cin, off_file_start_token);

    std::cin >> num_vertices >> num_faces >> num_edges;

    std::vector<Kernel::Point_2> vertices;
    std::vector<Kernel::Point_2> face_vertices;

    vertices.reserve(num_vertices);

    for(int vertex_index = 0; vertex_index != num_vertices; ++vertex_index){
        double x = 0.0;
        double y = 0.0;

        std::cin >> x >> y;

        vertices.emplace_back(x, y);
    }

    //Una sola faccia nel poligono semplice
    unsigned int num_vertices_in_face = 0;
    std::cin >> num_vertices_in_face;
    face_vertices.reserve(num_vertices_in_face);

    for(int i = 0; i != num_vertices_in_face; ++i){
        unsigned int face_vertex_index = 0;
        std::cin >> face_vertex_index;
        face_vertices.push_back(vertices[face_vertex_index]);
    }

    //Identifichiamo il vertice di partenza p0
    Point_2 p0 = *face_vertices.begin();
    Vertex_handle vh0 = dcel.insert_in_face_interior(p0, dcel.unbounded_face());

    //Il vertice p1 sarà l'ultimo vertice aggiunto alla dcel
    Point_2 p1 = p0;
    Vertex_handle vh1 = vh0;

    for(int i = 1; i < face_vertices.size(); ++i){
        Point_2 p = face_vertices[i];
        Vertex_handle vh = dcel.insert_in_face_interior(p, dcel.unbounded_face());

        dcel.insert_at_vertices(Segment_2(p1, p), vh1, vh);

        p1 = p;
        vh1 = vh;
    }

    dcel.insert_at_vertices(Segment_2(p0, p1), vh1, vh0);

    return dcel;
}
void write_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    const std::string off_file_start_token = "OFF";
    std::cout << off_file_start_token << "\n";

    unsigned int num_vertices = dcel.number_of_vertices();
    unsigned int num_faces = dcel.number_of_faces();
    unsigned int num_edges = dcel.number_of_edges();

    std::cout << num_vertices << " " << num_faces << " " << num_edges << "\n";

    //Mappiamo ogni indice al suo vertice
    std::map<Point_2, int> vertices;
    int vertex_index = 0;
    for(const auto& vh : dcel.vertex_handles()){
        vertices[vh->point()] = vertex_index;
        std::cout << vh->point().x() << " " << vh->point().y() << "\n";
        vertex_index += 1;
    }

    for(const auto& face : dcel.face_handles()){
        if(face->is_unbounded()) continue;

        int n_of_edges_for_face = 0;
        auto edge_circulator = face->outer_ccb();
        do {
            n_of_edges_for_face += 1;
        } while (++edge_circulator != face->outer_ccb());

        std::cout << n_of_edges_for_face << " ";

        edge_circulator = face->outer_ccb();
        do {
            std::cout << vertices[edge_circulator->source()->point()] << " ";
        } while (++edge_circulator != face->outer_ccb());

        std::cout << "\n";
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
        std::cerr << "Insert diagonal ";
        std::cerr << u->point().x() << " " << u->point().y() << " -> ";
        std::cerr << v->point().x() << " " << v->point().y() << "";
        std::cerr << "\n";
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

    const auto& u = helpers.at(b);
    if(types.at(u) == event_type::merge || types.at(v) == event_type::split){
        //Aggiunta della diagonale uv
        const Segment_2 segment(u->point(), v->point());
        std::cerr << "Insert diagonal ";
        std::cerr << u->point().x() << " " << u->point().y() << " -> ";
        std::cerr << v->point().x() << " " << v->point().y() << "";
        std::cerr << "\n";
        dcel.insert_at_vertices(segment, u, v);
    }

    helpers[b] = v;
}

//Procedura che identifica il tipo di evento relativo al vertice target di edge_handle
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

void print_edge(Halfedge_handle edge){
    std::cerr << edge->source()->point().x() << " " << edge->source()->point().y() << " -> ";
    std::cerr << edge->target()->point().x() << " " << edge->target()->point().y() << "  \n";
}

void print_vertex(Vertex_handle v){
    std::cerr << v->point().x() << " " << v->point().y() << "  \n";
}

void monotone_subdivision(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    //Per ogni vertice nella DCEL registriamo il tipo di evento
    std::map<Vertex_const_handle, event_type> types;

    std::vector<vertex_event> vertex_events;
    vertex_events.reserve(dcel.number_of_vertices());

    std::cout << dcel.number_of_faces() << "\n";
    Face_handle f = (++dcel.faces_begin());

    std::vector<Halfedge_handle> outer_edges;
    auto current = f->outer_ccb();
    do {
        Halfedge_handle he = current;
        outer_edges.push_back(he);
    }while (++current != f->outer_ccb());

    std::vector<Halfedge_handle> inner_edges;

    for(auto hi = f->holes_begin(); hi != f->holes_end(); ++hi){
        auto current_inner = *hi;
        do {
            Halfedge_handle he = current_inner;
            inner_edges.push_back(he);
        }while (++current_inner != *hi);
    }

    //Identificazione degli eventi associati a ciascun vertice
    //Per processare tutti i vertici una volta sola processiamo un edge alla volta
    //tenendo come vertice la sorgente dell'edge
    //Gli edge hanno puntatori ai vertici ma non viceversa
    for(auto& edge_handle : outer_edges){
        auto handle = edge_handle->twin();
        auto edge = *(handle);

        auto vertex = edge.source();
        const auto next_vertex = edge.target();
        const auto prev_vertex = edge.prev()->source();
        const auto point = vertex->point();
        const auto next_point = next_vertex->point();
        const auto prev_point = prev_vertex->point();

        const auto event_type = identify_event(prev_point, point, next_point);

        vertex_events.emplace_back(
                vertex,
                handle->prev()->twin(),
                handle->twin()
        );
        std::cerr << "Outer edge \n";
        print_vertex(vertex);
        print_edge(handle->prev()->twin());
        print_edge(handle->twin());

        types[vertex] = event_type;
    }

    for(auto& edge_handle : inner_edges){
        auto inverse_handle = edge_handle;

        auto edge = *(inverse_handle);

        auto vertex = edge.source();
        const auto next_vertex = edge.target();
        const auto prev_vertex = edge.prev()->source();

        const auto point = vertex->point();
        const auto next_point = next_vertex->point();
        const auto prev_point = prev_vertex->point();

        const auto event_type = identify_event(next_point, point, prev_point);

        vertex_events.emplace_back(
                vertex,
                inverse_handle,
                inverse_handle->prev()
        );
        std::cerr << "Inner edges \n";
        print_vertex(vertex);
        print_edge(inverse_handle);
        print_edge(inverse_handle->prev());

        types[vertex] = event_type;
    }

    //Ordiniamo gli eventi da sinistra a destra
    std::sort(vertex_events.begin(), vertex_events.end(), [](const auto& a, const auto& b){
        return CGAL::lexicographically_xy_smaller(a.vertex_handle->point(), b.vertex_handle->point());
    });

    std::for_each(vertex_events.begin(), vertex_events.end(), [](const vertex_event& e){
        std::cout << e.vertex_handle->point().x() << " " << e.vertex_handle->point().y() << "\n";
    });

    std::map<Halfedge_const_handle, Vertex_handle> helpers;
    EdgeStructure sweep_line_edges;

    //Processiamo ciascun evento
    for(const auto& event : vertex_events) {

        const auto& event_type = types[event.vertex_handle];

        if (event_type == event_type::start) {
            std::cerr << "event_type::start\n";
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::end) {
            std::cerr << "event_type::end\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
        }
        if (event_type == event_type::lower_regular) {
            std::cerr << "event_type::lower_regular\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::upper_regular) {
            std::cerr << "event_type::upper_regular\n";
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }
        if (event_type == event_type::split) {
            std::cerr << "event_type::split\n";
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
            std::cerr << "Insert next edge ";
            print_edge(event.next_edge);
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event_type == event_type::merge) {
            std::cerr << "event_type::merge\n";
            std::cerr << "Remove prev edge ";
            print_edge(event.prev_edge);
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            std::cerr << "Process next edge ";
            print_edge(event.next_edge);
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }

    }

}

int main() {

    //auto dcel = read_dcel();

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

    //print_dcel(dcel);
    monotone_subdivision(dcel);
    print_dcel(dcel);

    //write_dcel(dcel);

    return 0;
}