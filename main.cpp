#include <iostream>
#include <fstream>
#include <queue>

#include "OFFReader.cpp"

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
    event_type event_type;
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
    Halfedge_handle find_edge_underneath(const Halfedge_handle& e){
        sl = Line_2(e->source()->point(), Direction_2(0, 1));
        const auto it = edges.upper_bound(e);
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
        //Add diagonal to dcel from u to v
        std::cout << "Inserting edge: " << u->point().x() << " " << u->point().y() << " -> " << v->point().x() << " " << v->point().y() << "\n";
        Segment_2 segment(u->point(), v->point());
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
        //Add diagonal
        std::cout << "Inserting edge: " << u->point().x() << " " << u->point().y() << " -> " << v->point().x() << " " << v->point().y() << "\n";
        Segment_2 segment(u->point(), v->point());
        dcel.insert_at_vertices(segment, u, v);
    }

    helpers[b] = v;
}

vertex_event identify_event(const Halfedge_handle& edge_handle){
    auto edge = *(edge_handle->twin().ptr());

    auto vertex = edge.source();
    const auto next_vertex = edge.target();
    const auto prev_vertex = edge.prev()->source();

    const auto point = vertex->point();
    const auto next_point = next_vertex->point();
    const auto prev_point = prev_vertex->point();

    vertex_event event;
    event.vertex_handle = vertex;
    event.next_edge = edge_handle->twin();
    event.prev_edge = edge_handle->twin()->prev();

    if (CGAL::lexicographically_xy_smaller(prev_point, point) &&
        CGAL::lexicographically_xy_smaller(point, next_point)) {
        //lower regular
        event.event_type = event_type::lower_regular;
    }
    if (CGAL::lexicographically_xy_larger(prev_point, point) &&
        CGAL::lexicographically_xy_larger(point, next_point)) {
        //upper regular
        event.event_type = event_type::upper_regular;
    }
    if (CGAL::lexicographically_xy_larger(prev_point, point) &&
        CGAL::lexicographically_xy_larger(next_point, point)) {
        //Based on turn direction identify start or split event
        if (CGAL::left_turn(prev_point, point, next_point)) {
            event.event_type = event_type::start;
        } else {
            event.event_type = event_type::split;
        }
    }
    if (CGAL::lexicographically_xy_smaller(prev_point, point) &&
        CGAL::lexicographically_xy_smaller(next_point, point)) {
        //Based on turn direction identify end or merge event
        if (CGAL::left_turn(prev_point, point, next_point)) {
            event.event_type = event_type::end;
        } else {
            event.event_type = event_type::merge;
        }
    }

    return event;
}

void process_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    //Per ogni vertice nella DCEL registriamo il tipo di evento
    std::map<Vertex_const_handle, event_type> types;

    std::vector<vertex_event> vertex_events;
    vertex_events.reserve(dcel.number_of_vertices());

    //Identificazione degli eventi associati a ciascun vertice
    //Per processare tutti i vertici una volta sola processiamo un edge alla volta
    //tenendo come vertice la sorgente dell'edge
    //Gli edge hanno puntatori ai vertici ma non viceversa
    for(auto& edge_handle : dcel.edge_handles()){
        const auto& event = identify_event(edge_handle);

        vertex_events.push_back(event);
        types[edge_handle->source()] = event.event_type;
    }

    std::sort(vertex_events.begin(), vertex_events.end(), [](const auto& a, const auto& b){
        return CGAL::lexicographically_xy_smaller(a.vertex_handle->point(), b.vertex_handle->point());
    });

    std::map<Halfedge_const_handle, Vertex_handle> helpers;

    EdgeStructure sweep_line_edges;

    for(const auto& event : vertex_events) {

        if (event.event_type == event_type::start) {
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event.event_type == event_type::end) {
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
        }
        if (event.event_type == event_type::lower_regular) {
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event.event_type == event_type::upper_regular) {
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }
        if (event.event_type == event_type::split) {
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event.event_type == event_type::merge) {
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            process(event.vertex_handle, event.next_edge, sweep_line_edges, helpers, types, dcel);
        }

    }

}

void print_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    std::cout << "Printing DCEL\n";

    std::cout << "Vertices\n";
    for(Vertex_handle& v : dcel.vertex_handles()){
        std::cout << v->point().x() << " " << v->point().y() << "\n";
    }

    std::cout << "Half-edges\n";
    for(const auto& he : dcel.edge_handles()){
        std::cout << he->twin()->source()->point().x() << " " << he->twin()->source()->point().y() << " -> " << he->twin()->target()->point().x() << " " << he->twin()->target()->point().y() << "\n";
    }

}

void output_dcel(const std::string_view filename, CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    std::ofstream off_file(filename.data());

    const std::string off_file_start_token = "OFF";
    off_file << off_file_start_token << "\n";

    unsigned int num_vertices = dcel.number_of_vertices();
    unsigned int num_faces = dcel.number_of_faces();
    unsigned int num_edges = dcel.number_of_edges();

    off_file << num_vertices << " " << num_faces << " " << num_edges << "\n";

    //Map each vertex to its index
    std::map<Point_2, int> vertices;
    int vertex_index = 1;
    for(const auto& vh : dcel.vertex_handles()){
        vertices[vh->point()] = vertex_index;
        off_file << vh->point().x() << " " << vh->point().y() << "\n";
        vertex_index += 1;
    }

    for(const auto& face : dcel.face_handles()){
        if(face->is_unbounded()) continue;

        int n_of_edges_for_face = 0;
        auto edge_circulator = face->outer_ccb();
        do {
            n_of_edges_for_face += 1;
        } while (++edge_circulator != face->outer_ccb());

        off_file << n_of_edges_for_face << " ";

        edge_circulator = face->outer_ccb();
        do {
            off_file << vertices[edge_circulator->source()->point()] << " ";
        } while (++edge_circulator != face->outer_ccb());

        off_file << "\n";
    }

}

int main() {

    //Example simple polygon
    Geometry polygon = OFFReader::read("cube.off");

    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    //Creazione dei segmenti
    std::vector<Segment_2> dcel_segments;
    dcel_segments.reserve(polygon.edges.size());
    std::transform(polygon.edges.begin(), polygon.edges.end(), std::back_inserter(dcel_segments), [&](const auto& pair){
        const auto a = polygon.vertices[pair.first];
        const auto b = polygon.vertices[pair.second];
        return Segment_2(a, b);
    });

    //Posizionamento di tutti i vertici nella faccia infinita
    std::map<Point_2, Vertex_handle> dcel_vertices;
    for(const auto& vertex : polygon.vertices){
        //Salviamo i vertex_handle in una mappa per connetterli facilmente con i segmenti
        dcel_vertices[vertex] = dcel.insert_in_face_interior(vertex, dcel.unbounded_face());
    }
    for(auto & dcel_segment : dcel_segments){
        //Per ogni segmento connettiamo i corrispondenti vertici
        dcel.insert_at_vertices(dcel_segment, dcel_vertices[dcel_segment.source()], dcel_vertices[dcel_segment.target()]);
    }

    print_dcel(dcel);
    std::cout << "--------------------- \n";
    process_dcel(dcel);

    print_dcel(dcel);

    output_dcel("output.off", dcel);

    return 0;
}