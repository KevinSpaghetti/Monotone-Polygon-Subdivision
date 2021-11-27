#include <iostream>
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

struct sweep_line_event {
    Vertex_handle vertex_handle;
    Halfedge_handle next_edge;
    Halfedge_handle prev_edge;
    event_type event_type;
};

class EdgeStructure {
private:
    std::set<Halfedge_handle> edges;

public:

    void insert(const Halfedge_handle& e){
        edges.insert(e);
    }
    void erase(const Halfedge_handle& e){
        edges.erase(e);
    }

    Halfedge_handle find_edge_underneath(const Vertex_handle& v){

        Halfedge_handle e;
        float min_distance = std::numeric_limits<float>::max();

        for(const auto edge : edges){

            const auto line = Line_2(v->point(), Direction_2(0, 1));
            CGAL::Object u = CGAL::intersection(edge->curve().line(), line);
            if(CGAL::object_cast<Point_2>(&u)){
                Point_2 intersection_point = *CGAL::object_cast<Point_2>(&u);
                if(intersection_point.y() < v->point().y()){
                    if(v->point().y() - intersection_point.y() < min_distance){
                        min_distance = v->point().y() - intersection_point.y();
                        e = edge;
                    }
                }
            }

        }

        return e;
    }

    void print_content(){
        for(const auto e : edges){
            std::cout << e->source()->point().x() << " " << e->source()->point().y();
            std::cout << " -> ";
            std::cout << e->target()->point().x() << " " << e->target()->point().y() << "\n";
            std::cout << std::endl;
        }
    }

};


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
        Segment_2 segment(u->point(), v->point());
        dcel.insert_at_vertices(segment, u, v);
    }

    sweep_line_edges.erase(e);
}
void process(const Vertex_handle& v,
             EdgeStructure& sweep_line_edges,
             std::map<Halfedge_const_handle, Vertex_handle>& helpers,
             const std::map<Vertex_const_handle, event_type>& types,
             CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    const auto b = sweep_line_edges.find_edge_underneath(v);

    const auto& u = helpers.at(b);
    if(types.at(u) == event_type::merge || types.at(v) == event_type::split){
        //Add diagonal
        Segment_2 segment(u->point(), v->point());
        dcel.insert_at_vertices(segment, u, v);
    }

    helpers[b] = v;
}

void process_dcel(CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>>& dcel){

    //Monotone subdivision of polygon
    std::map<Vertex_const_handle, event_type> types;
    //Events
    const auto sweep_line_event_comparator = [](const auto &a, const auto &b) {
        return CGAL::lexicographically_xy_larger(a.vertex_handle->point(), b.vertex_handle->point());
    };
    std::priority_queue<sweep_line_event, std::vector<sweep_line_event>, decltype(sweep_line_event_comparator)> vertex_events(
            sweep_line_event_comparator);

    for(auto& edge_handle : dcel.edge_handles()){

        auto edge = *(edge_handle->twin().ptr());

        auto vertex = edge.source();
        const auto next_vertex = edge.target();
        const auto prev_vertex = edge.prev()->source();

        const auto point = vertex->point();
        const auto next_point = next_vertex->point();
        const auto prev_point = prev_vertex->point();

        sweep_line_event event;
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

        vertex_events.push(event);
        types[vertex] = event.event_type;
    }

    std::map<Halfedge_const_handle, Vertex_handle> helpers;

    Line_2 sl;
    EdgeStructure sweep_line_edges;

    while (!vertex_events.empty()) {
        auto &event = vertex_events.top();

        const auto prev_edge_source_point = event.prev_edge->source()->point();
        const auto prev_edge_target_point = event.prev_edge->target()->point();

        const auto next_edge_source_point = event.next_edge->source()->point();
        const auto next_edge_target_point = event.next_edge->target()->point();

        sweep_line_edges.print_content();

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
            process(event.vertex_handle, sweep_line_edges, helpers, types, dcel);
        }
        if (event.event_type == event_type::split) {
            process(event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            insert(event.next_edge, sweep_line_edges, helpers, types);
        }
        if (event.event_type == event_type::merge) {
            remove(event.prev_edge, event.vertex_handle, sweep_line_edges, helpers, types, dcel);
            process(event.vertex_handle, sweep_line_edges, helpers, types, dcel);
        }

        vertex_events.pop();
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

int main() {

    //Example simple polygon
    Geometry polygon = OFFReader::read("cube.off");

    CGAL::Arrangement_2<CGAL::Arr_segment_traits_2<Kernel>> dcel;

    std::vector<Segment_2> dcel_segments;

    dcel_segments.reserve(polygon.edges.size());

    std::transform(polygon.edges.begin(), polygon.edges.end(), std::back_inserter(dcel_segments), [&](const auto& pair){
        const auto a = polygon.vertices[pair.first];
        const auto b = polygon.vertices[pair.second];
        return Segment_2(a, b);
    });

    std::map<Point_2, Vertex_handle> dcel_vertices;
    for(const auto& vertex : polygon.vertices){
        dcel_vertices[vertex] = dcel.insert_in_face_interior(vertex, dcel.unbounded_face());
    }
    for(int i = 0; i < dcel_segments.size(); ++i){
        dcel.insert_at_vertices(dcel_segments[i], dcel_vertices[dcel_segments[i].source()], dcel_vertices[dcel_segments[i].target()]);
    }

    print_dcel(dcel);
    std::cout << "--------------------- \n";
    process_dcel(dcel);

    print_dcel(dcel);

    return 0;
}