#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <set>
#include <fstream>
#include <CGAL/number_utils.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Plane_3 Plane;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Color {
    float r, g, b;
};

// mouse rotation
float rotation_x = 0.0f;
float rotation_y = 0.0f;
double last_cursor_x, last_cursor_y;
bool mouse_pressed = false;

bool is_degenerate_tetrahedron(const Point& p0, const Point& p1, const Point& p2, const Point& p3) {
    return CGAL::coplanar(p0, p1, p2, p3);
}

Nef_polyhedron make_tetrahedron(Point center, double size) {
    double s = size;
    std::vector<Point> points = {
        Point(center.x()+s, center.y()+s, center.z()+s),
        Point(center.x()-s, center.y()-s, center.z()+s),
        Point(center.x()-s, center.y()+s, center.z()-s),
        Point(center.x()+s, center.y()-s, center.z()-s),
    };
    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(4);
    for (int i = 0; i < 4; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c) {
        mesh.add_face(vi[a], vi[b], vi[c]);
    };

    add_face(0, 2, 1);
    add_face(0, 1, 3);
    add_face(0, 3, 2);
    add_face(1, 2, 3);

    return Nef_polyhedron(mesh);
}

Nef_polyhedron make_cube(Point center, double size) {
    double s = size;
    std::vector<Point> points = {
        Point(center.x()-s, center.y()-s, center.z()-s),
        Point(center.x()+s, center.y()-s, center.z()-s),
        Point(center.x()+s, center.y()+s, center.z()-s),
        Point(center.x()-s, center.y()+s, center.z()-s),
        Point(center.x()-s, center.y()-s, center.z()+s),
        Point(center.x()+s, center.y()-s, center.z()+s),
        Point(center.x()+s, center.y()+s, center.z()+s),
        Point(center.x()-s, center.y()+s, center.z()+s),
    };
    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(8);
    for (int i = 0; i < 8; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c, int d) {
        mesh.add_face(vi[a], vi[b], vi[c]);
        mesh.add_face(vi[a], vi[c], vi[d]);
    };

    add_face(0, 1, 2, 3);
    add_face(7, 6, 5, 4);
    add_face(0, 4, 5, 1);
    add_face(1, 5, 6, 2);
    add_face(2, 6, 7, 3);
    add_face(3, 7, 4, 0);

    return Nef_polyhedron(mesh);
}


Nef_polyhedron make_Octahedron(Point center, double size) {
    double s = size*2;

    std::vector<Point> points = {
        Point(center.x(),     center.y(),     center.z() + s),
        Point(center.x(),     center.y(),     center.z() - s),
        Point(center.x() + s, center.y(),     center.z()),
        Point(center.x() - s, center.y(),     center.z()),
        Point(center.x(),     center.y() + s, center.z()),
        Point(center.x(),     center.y() - s, center.z())
    };

    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(6);
    for (int i = 0; i < 6; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c) {
        mesh.add_face(vi[a], vi[b], vi[c]);
    };

    add_face(0, 2, 4);
    add_face(0, 4, 3);
    add_face(0, 3, 5);
    add_face(0, 5, 2);
    add_face(1, 4, 2);
    add_face(1, 3, 4);
    add_face(1, 5, 3);
    add_face(1, 2, 5);

    return Nef_polyhedron(mesh);
}

Nef_polyhedron make_Dodecahedron(Point center, double size){
    double s = size;
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
     std::vector<Point> points = {
        Point(center.x()+s, center.y()+s, center.z()+s),
        Point(center.x()+s, center.y()+1.0*s, center.z()-1.0*s),
        Point(center.x()+s, center.y()-1.0*s, center.z()+1.0*s),
        Point(center.x()+s, center.y()-1.0*s, center.z()-1.0*s),
        Point(center.x()-1.0*s, center.y()+1.0*s, center.z()+1.0*s),
        Point(center.x()-1.0*s, center.y()+1.0*s, center.z()-1.0*s),
        Point(center.x()-1.0*s, center.y()-1.0*s, center.z()+1.0*s),
        Point(center.x()-1.0*s, center.y()-1.0*s, center.z()-1.0*s),
        Point(center.x(), center.y()+s/phi, center.z()+phi*s),
        Point(center.x(), center.y()+s/phi, center.z()-phi*s),
        Point(center.x(), center.y()-1.0*s/phi, center.z()+phi*s),
        Point(center.x(), center.y()-1.0*s/phi, center.z()-phi*s),
        Point(center.x()+s/phi, center.y()+phi*s, center.z()),
        Point(center.x()+s/phi, center.y()-phi*s, center.z()),
        Point(center.x()-1.0*s/phi, center.y()+phi*s, center.z()),
        Point(center.x()-1.0*s/phi, center.y()-phi*s, center.z()),
        Point(center.x()+phi*s, center.y(), center.z()+s/phi),
        Point(center.x()+phi*s, center.y(), center.z()-1.0*s/phi),
        Point(center.x()-phi*s, center.y(), center.z()+s/phi),
        Point(center.x()-phi*s, center.y(), center.z()-1.0*s/phi),
    };
    Surface_mesh mesh;
    std::vector<Surface_mesh::Vertex_index> vi(20);
    for (int i = 0; i < 20; ++i) vi[i] = mesh.add_vertex(points[i]);

    auto add_triangle_fan = [&](std::array<int, 5> f) {
        mesh.add_face(vi[f[0]], vi[f[1]], vi[f[2]]);
        mesh.add_face(vi[f[0]], vi[f[2]], vi[f[3]]);
        mesh.add_face(vi[f[0]], vi[f[3]], vi[f[4]]);
    };

    // 12 pentagonal faces divided into 3 triangles each
    std::vector<std::array<int,5>> faces = {
        {0, 8, 4, 14, 12},
        {0, 12, 1, 17, 16},
        {0, 16, 2, 10, 8},
        {8, 10, 6, 18, 4},
        {4, 18, 19, 5, 14},
        {14, 5, 9, 1, 12},
        {1, 9, 11, 3, 17},
        {17, 3, 13, 2, 16},
        {2, 13, 15, 6, 10},
        {6, 15, 7, 19, 18},
        {5, 19, 7, 11, 9},
        {3, 11, 7, 15, 13}
    };



    for (const auto& face : faces) {
        add_triangle_fan(face);
    }

    return Nef_polyhedron(mesh);
}

Nef_polyhedron make_Icosahedron(Point center, double size){
    double s = size;
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<Point> points = {
        Point(center.x(), center.y()+s, center.z()+phi*s),
        Point(center.x(), center.y()+s, center.z()-phi*s),
        Point(center.x(), center.y()-1.0*s, center.z()+phi*s),
        Point(center.x(), center.y()-1.0*s, center.z()-phi*s),
        Point(center.x()+s, center.y()+phi*s, center.z()),
        Point(center.x()+s, center.y()-phi*s, center.z()),
        Point(center.x()-1.0*s, center.y()+phi*s, center.z()),
        Point(center.x()-1.0*s, center.y()-phi*s, center.z()),
        Point(center.x()+phi*s, center.y(), center.z()+s),
        Point(center.x()+phi*s, center.y(), center.z()-1.0*s),
        Point(center.x()-phi*s, center.y(), center.z()+s),
        Point(center.x()-phi*s, center.y(), center.z()-1.0*s),
    };
    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(12);
    for (int i = 0; i < 12; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c) {
        mesh.add_face(vi[a], vi[b], vi[c]);
    };

    add_face(0, 4, 6);
    add_face(0, 6, 10);
    add_face(0, 10, 2);
    add_face(0, 2, 8);
    add_face(0, 8, 4);
    add_face(4, 8, 9);
    add_face(4, 9, 1);
    add_face(4, 1, 6);
    add_face(6, 1, 11);
    add_face(6, 11, 10);
    add_face(10, 11, 7);
    add_face(10, 7, 2);
    add_face(2, 7, 5);
    add_face(2, 5, 8);
    add_face(8, 5, 9);
    add_face(9, 5, 3);
    add_face(9, 3, 1);
    add_face(1, 3, 11);
    add_face(11, 3, 7);
    add_face(7, 3, 5);

    return Nef_polyhedron(mesh);
}



Surface_mesh rotate_mesh(const Surface_mesh& input_mesh, const Kernel::Vector_3& axis, double angle_degrees) {
    Surface_mesh output_mesh = input_mesh;
    Kernel::Vector_3 n = axis / std::sqrt(CGAL::to_double(axis.squared_length()));
    double angle_radians = angle_degrees * CGAL_PI / 180.0;
    double c = std::cos(angle_radians);
    double s = std::sin(angle_radians);
    double t = 1.0 - c;
    double x = CGAL::to_double(n.x()), y = CGAL::to_double(n.y()), z = CGAL::to_double(n.z());

    // Rodrigues' rotation formula
    CGAL::Aff_transformation_3<Kernel> rotation(
        x * x * t + c,     x * y * t - z * s, x * z * t + y * s,
        y * x * t + z * s, y * y * t + c,     y * z * t - x * s,
        z * x * t - y * s, z * y * t + x * s, z * z * t + c
    );
    for (auto v : output_mesh.vertices()) {
        output_mesh.point(v) = rotation.transform(output_mesh.point(v));
    }
    return output_mesh;
}

Nef_polyhedron rotate_nef_polyhedron(const Nef_polyhedron& nef, const Kernel::Vector_3& axis, double angle_degrees) {
    Surface_mesh mesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh);
    Surface_mesh rotated_mesh = rotate_mesh(mesh, axis, angle_degrees);
    return Nef_polyhedron(rotated_mesh);
}

void tag_original_edges(Surface_mesh& mesh) {
    auto edge_map_pair = mesh.add_property_map<Surface_mesh::Edge_index, bool>("e:original", false);
    auto& edge_map = edge_map_pair.first;
    for (auto edge : mesh.edges())
        edge_map[edge] = true;
}


void set_perspective(float fovy, float aspect, float zNear, float zFar) {
    float f = 1.0f / tanf(fovy * 0.5f * 3.14159265f / 180.0f);
    float m[16] = { 0 };
    m[0] = f / aspect;
    m[5] = f;
    m[10] = (zFar + zNear) / (zNear - zFar);
    m[11] = -1.0f;
    m[14] = (2.0f * zFar * zNear) / (zNear - zFar);
    glMultMatrixf(m);
}

void assign_random_face_colors(Surface_mesh& mesh) {
    using Face_index = Surface_mesh::Face_index;
    auto color_map_pair = mesh.add_property_map<Face_index, Color>("f:color", Color(0, 0, 0));
    auto& color_map = color_map_pair.first;
    for (auto face : mesh.faces())
        color_map[face] = Color(
            static_cast<unsigned char>(rand() % 256),
            static_cast<unsigned char>(rand() % 256),
            static_cast<unsigned char>(rand() % 256)
        );
}


void triangulate_faces_with_colors(Surface_mesh& mesh) {
    tag_original_edges(mesh);
    assign_random_face_colors(mesh);
    using Face_index = Surface_mesh::Face_index;
    using Vertex_index = Surface_mesh::Vertex_index;
    auto color_map_opt = mesh.property_map<Face_index, Color>("f:color");
    if (!color_map_opt.has_value())
        throw std::runtime_error("Missing 'f:color' map.");
    auto& color_map = *color_map_opt;
    std::map<Face_index, Color> face_colors;
    std::vector<std::tuple<Point, Point, Point, Color>> triangles;
    for (auto face : mesh.faces()) {
        Color col = color_map[face];
        std::vector<Point> points;
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(face), mesh))
            points.push_back(mesh.point(v));
        if (points.size() >= 3)
            for (std::size_t i = 1; i + 1 < points.size(); ++i)
                triangles.emplace_back(points[0], points[i], points[i + 1], col);
    }
    std::set<std::pair<Point, Point>> original_edges;
    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));
        if (p1 > p2) std::swap(p1, p2); // Ensure consistent ordering
        original_edges.insert({p1, p2});
    }
    mesh.clear();
    auto new_color_map_pair = mesh.add_property_map<Face_index, Color>("f:color", Color(0, 0, 0));
    auto& new_color_map = new_color_map_pair.first;
    std::map<Point, Vertex_index> point_map;
    auto add_vertex = [&](const Point& p) -> Vertex_index {
        auto it = point_map.find(p);
        if (it != point_map.end()) return it->second;
        Vertex_index v = mesh.add_vertex(p);
        point_map[p] = v;
        return v;
    };
    for (const auto& [p0, p1, p2, col] : triangles) {
        Vertex_index v0 = add_vertex(p0);
        Vertex_index v1 = add_vertex(p1);
        Vertex_index v2 = add_vertex(p2);
        Face_index f = mesh.add_face(v0, v1, v2);
        if (f != Surface_mesh::null_face()) new_color_map[f] = col;
    }
    auto new_edge_map_pair = mesh.add_property_map<Surface_mesh::Edge_index, bool>("e:original", false);
    auto& new_edge_map = new_edge_map_pair.first;
    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));
        if (p1 > p2) std::swap(p1, p2);
        if (original_edges.count({p1, p2}) > 0)
            new_edge_map[edge] = true;
    }
}


void draw_mesh(const Surface_mesh& mesh,
    const std::map<Surface_mesh::Face_index, Color>& face_colors,
    const std::map<Surface_mesh::Edge_index, Color>& edge_colors,
    const std::set<Surface_mesh::Edge_index>& reflex_edges_set)
    {
    // draw faces
    for (auto face : mesh.faces()) {
        glBegin(GL_TRIANGLE_FAN);
        auto color = face_colors.at(face);
        glColor3f(color.r, color.g, color.b);

        for (auto v : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            Point p = mesh.point(v);
            glVertex3f(CGAL::to_double(p.x()),
                       CGAL::to_double(p.y()),
                       CGAL::to_double(p.z()));
        }
        glEnd();
    }

    // draw edges
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for (auto edge : mesh.edges()) {
        auto color = edge_colors.at(edge);   // get color
        glColor3f(color.r, color.g, color.b);

        bool isReflexEdge = reflex_edges_set.find(edge) != reflex_edges_set.end();

        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));

        // reflex edge only slightly toward the front
        if (isReflexEdge) {
            auto f1 = mesh.face(h);
            auto f2 = mesh.face(mesh.opposite(h));

            Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
            Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);
            Kernel::Vector_3 avg_n = (n1 + n2) / 2.0;

            double len = std::sqrt(CGAL::to_double(avg_n.squared_length()));
            if (len > 1e-8) avg_n = avg_n / len;

            double offset = 0.1;  // lift slightly
            p1 = Point(p1.x() + avg_n.x() * offset,
                       p1.y() + avg_n.y() * offset,
                       p1.z() + avg_n.z() * offset);
            p2 = Point(p2.x() + avg_n.x() * offset,
                       p2.y() + avg_n.y() * offset,
                       p2.z() + avg_n.z() * offset);
        }

        glLineWidth(isReflexEdge ? 8.0f : 3.0f); // reflex is thick
        glVertex3f(CGAL::to_double(p1.x()),
                   CGAL::to_double(p1.y()),
                   CGAL::to_double(p1.z()));
        glVertex3f(CGAL::to_double(p2.x()),
                   CGAL::to_double(p2.y()),
                   CGAL::to_double(p2.z()));
    }
    glEnd();
}

void draw_skeleton(const Surface_mesh& mesh,
                   float vertex_size = 4.0,
                   float edge_width = 2.0)
{
    // edges
    glLineWidth(edge_width);
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));

        glVertex3f(CGAL::to_double(p1.x()),
                   CGAL::to_double(p1.y()),
                   CGAL::to_double(p1.z()));
        glVertex3f(CGAL::to_double(p2.x()),
                   CGAL::to_double(p2.y()),
                   CGAL::to_double(p2.z()));
    }
    glEnd();

    // vertices
    glPointSize(vertex_size);
    glBegin(GL_POINTS);
    glColor3f(1.0f, 0.3f, 0.3f); // red
    for (auto v : mesh.vertices()) {
        Point p = mesh.point(v);
        glVertex3f(CGAL::to_double(p.x()),
                   CGAL::to_double(p.y()),
                   CGAL::to_double(p.z()));
    }
    glEnd();
}

void draw_guard_edges(const Surface_mesh& mesh,
                      const std::set<Surface_mesh::Edge_index>& guard_edges)
{
    glLineWidth(5.0f);
    glColor3f(1.0f, 1.0f, 0.0f); // yellow
    glBegin(GL_LINES);
    for (auto edge : guard_edges) {
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));

        glVertex3f(CGAL::to_double(p1.x()),
                   CGAL::to_double(p1.y()),
                   CGAL::to_double(p1.z()));
        glVertex3f(CGAL::to_double(p2.x()),
                   CGAL::to_double(p2.y()),
                   CGAL::to_double(p2.z()));
    }
    glEnd();
}

// mouse callback
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    if (mouse_pressed) {
        float dx = static_cast<float>(xpos - last_cursor_x);
        float dy = static_cast<float>(ypos - last_cursor_y);
        rotation_y += dx * 0.2f;
        rotation_x += dy * 0.2f;
    }
    last_cursor_x = xpos;
    last_cursor_y = ypos;
}

// mouse button callback
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS)
            mouse_pressed = true;
        else if (action == GLFW_RELEASE)
            mouse_pressed = false;
    }
}


//1-skelton
std::map<Surface_mesh::Vertex_index, std::set<Surface_mesh::Vertex_index>> One_skelton(const Surface_mesh& mesh) {
    std::map<Surface_mesh::Vertex_index, std::set<Surface_mesh::Vertex_index>> adjacency;

    for (auto edge : mesh.edges()) {
        auto h = mesh.halfedge(edge);
        auto v1 = mesh.source(h);
        auto v2 = mesh.target(h);
        adjacency[v1].insert(v2);
        adjacency[v2].insert(v1);
    }

    return adjacency;
}



//Greedy algorithm

std::set<Surface_mesh::Edge_index> GreedyEdgeGuards(
    const Surface_mesh& mesh,
    const std::set<Surface_mesh::Edge_index>& diagonal_edges_set)
{
    using V = Surface_mesh::Vertex_index;
    using E = Surface_mesh::Edge_index;

    std::set<V> S;      // Covered vertices
    std::set<E> L;      // Selected guard edges
    
    std::vector<V> Vlist;
    for (auto v : mesh.vertices())
        Vlist.push_back(v);

    // Build adjacency (1-ring neighbors) and incident edges
    std::map<V, std::vector<V>> adj;
    std::map<V, std::vector<E>> incident_edges;

    for (auto e : mesh.edges()) {
        auto h = mesh.halfedge(e);
        auto v1 = mesh.source(h);
        auto v2 = mesh.target(h);

        adj[v1].push_back(v2);
        adj[v2].push_back(v1);

        incident_edges[v1].push_back(e);
        incident_edges[v2].push_back(e);
    }

    // Greedy loop
    for (auto v : Vlist) {

        if (S.count(v)) continue;

        S.insert(v);
    
        E new_edge;
        bool added=false;
        for (auto e : incident_edges[v]) {
            auto h = mesh.halfedge(e);
            auto a = mesh.source(h);
            auto b = mesh.target(h);
            V u = (a == v ? b : a);
            if (diagonal_edges_set.count(e)) continue;
            new_edge=e;
            if (S.count(u)) continue;
            S.insert(u);
            L.insert(e);
            added=true;
            break;
        }
        if(!added)L.insert(new_edge);
    }

    return L;
}



bool is_same_polygon_face(const Surface_mesh& mesh, Surface_mesh::Face_index f1, Surface_mesh::Face_index f2) {
    if (f1 == Surface_mesh::null_face() || f2 == Surface_mesh::null_face())
        return false;

    Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
    Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

    // The normal vectors are almost parallel
    return std::fabs(CGAL::to_double(n1 * n2 - 1.0)) < 1e-6;
}


void classify_edges_ver3(const Surface_mesh& mesh,
                         std::map<Surface_mesh::Edge_index, Color>& edge_colors,
                         std::set<Surface_mesh::Edge_index>& reflex_edges_set,
                         std::set<Surface_mesh::Edge_index>& diagonal_edges_set) {
    int real_edges = 0;
    int reflex_edges = 0;
    int convex_edges = 0;
    int diagonal_edges = 0;

    int classA_convex = 0, classA_reflex = 0;
    int classB_convex = 0, classB_reflex = 0;
    int classC_convex = 0, classC_reflex = 0;
    int classD_convex = 0, classD_reflex = 0;

    for (auto edge : mesh.edges()) {
        auto h1 = mesh.halfedge(edge, 0);

        // skip border
        if (mesh.is_border(h1)) {
            real_edges++;
            edge_colors[edge] = {1.0f, 1.0f, 1.0f}; // white
            continue;
        }

        auto h2 = mesh.opposite(h1);
        auto f1 = mesh.face(h1);
        auto f2 = mesh.face(h2);

        // normal
        Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
        Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

        // check diagonal
        double dot = CGAL::to_double(n1 * n2);
        if (std::fabs(dot - 1.0) < 1e-6) {
            diagonal_edges++;
            edge_colors[edge] = {0.7f, 0.7f, 0.7f};
            diagonal_edges_set.insert(edge);
            continue;
        }

        // check reflex or convex
        auto v_source = mesh.source(h1);
        auto v_other = mesh.target(mesh.next(h2));
        Kernel::Vector_3 test_vec = mesh.point(v_other) - mesh.point(v_source);

        bool isReflex = (test_vec * n1 > 0);
        if (isReflex){
            reflex_edges++;
            reflex_edges_set.insert(edge);  // insert reflex edge
        }
        else
            convex_edges++;
        real_edges++;

        // classified
        double x1 = CGAL::to_double(n1.x());
        double y1 = CGAL::to_double(n1.y());
        double x2 = CGAL::to_double(n2.x());
        double y2 = CGAL::to_double(n2.y());

        bool isA = (y1 >= 0 && y2 >= 0);
        bool isB = (y1 < 0 && y2 < 0);
        auto vA  = mesh.source(h1);
        auto vB  = mesh.target(h1);
        auto vC1 = mesh.target(mesh.next(h1));
        auto vC2 = mesh.target(mesh.next(h2));
        Kernel::Point_3 A = mesh.point(vA);
        Kernel::Point_3 B = mesh.point(vB);
        Kernel::Point_3 C1 = mesh.point(vC1);
        Kernel::Point_3 C2 = mesh.point(vC2);
        Kernel::Vector_3 up(0.0, 1.0, 0.0);
        Kernel::Vector_3 e = B - A;
        Kernel::Vector_3 m = CGAL::cross_product(e, up);
        double s1 = CGAL::to_double(m * (C1 - A));
        double s2 = CGAL::to_double(m * (C2 - A));
        bool isC = (s1 < 0.0 && s2 < 0.0);
        bool isD = (s1 > 0.0 && s2 > 0.0);


        // count and color
        if (isA) {
            if (isReflex) { classA_reflex++; edge_colors[edge] = {1.0f, 0.5f, 0.5f}; } // Reflex A2(pale red)
            else          { classA_convex++; edge_colors[edge] = {1.0f, 0.0f, 0.0f}; } // Convex A1(red)
        }
        else if (isB) {
            if (isReflex) { classB_reflex++; edge_colors[edge] = {0.5f, 0.5f, 1.0f}; } // Reflex B2(pale blue)
            else          { classB_convex++; edge_colors[edge] = {0.0f, 0.0f, 1.0f}; } // Convex B1(blue)
        }
        else if (isC) {
            if (isReflex) { classC_reflex++; edge_colors[edge] = {0.5f, 1.0f, 0.5f}; } // Reflex C2(pale green)
            else          { classC_convex++; edge_colors[edge] = {0.0f, 1.0f, 0.0f}; } // Convex C1(green)
        }
        else if (isD) {
            if (isReflex) { classD_reflex++; edge_colors[edge] = {1.0f, 1.0f, 0.5f}; } // Reflex D2(pale yellow)
            else          { classD_convex++; edge_colors[edge] = {1.0f, 1.0f, 0.0f}; } // Convex D1(yellow)
        }
        else {
            edge_colors[edge] = {1.0f, 1.0f, 1.0f};
        }
    }

    // result
    std::cout << "[Edge Classification Result]\n";
    std::cout << "Total Edges    : " << mesh.number_of_edges() << "\n";
    std::cout << "Real Edges     : " << real_edges << "\n";
    std::cout << "Reflex Edges   : " << reflex_edges << "\n";
    std::cout << "Convex Edges   : " << convex_edges << "\n";
    std::cout << "Diagonal Edges : " << diagonal_edges << "\n\n";

    std::cout << "[A to D Reflex/Convex Classes]\n";
    std::cout << "A1 (Y>=0 both) Convex : " << classA_convex << "\n";
    std::cout << "A2 (Y>=0 both) Reflex : " << classA_reflex << "\n";
    std::cout << "B1 (Y<0 both) Convex : " << classB_convex << "\n";
    std::cout << "B2 (Y<0 both) Reflex : " << classB_reflex << "\n";
    std::cout << "C1 (X>=0 or mixed +Y/-Y) Convex : " << classC_convex << "\n";
    std::cout << "C2 (X>=0 or mixed +Y/-Y) Reflex : " << classC_reflex << "\n";
    std::cout << "D1 (X<0 or mixed +Y/-Y) Convex : " << classD_convex << "\n";
    std::cout << "D2 (X<0 or mixed +Y/-Y) Reflex : " << classD_reflex << "\n\n";
}

char classify_edge_class(
    const Surface_mesh& mesh,
    Surface_mesh::Edge_index e)
{
    auto h1 = mesh.halfedge(e, 0);

    // border is skip
    if (mesh.is_border(h1)) return 'X';

    auto h2 = mesh.opposite(h1);
    auto f1 = mesh.face(h1);
    auto f2 = mesh.face(h2);

    // normal vector
    Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
    Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

    double x1 = CGAL::to_double(n1.x());
    double y1 = CGAL::to_double(n1.y());
    double x2 = CGAL::to_double(n2.x());
    double y2 = CGAL::to_double(n2.y());
    

    // class A
    if (y1 >= 0 && y2 >= 0)
        return 'A';

    // class B
    if (y1 < 0 && y2 < 0)
        return 'B';
    auto vA  = mesh.source(h1);
    auto vB  = mesh.target(h1);
    auto vC1 = mesh.target(mesh.next(h1));
    auto vC2 = mesh.target(mesh.next(h2));
    Kernel::Point_3 A = mesh.point(vA);
    Kernel::Point_3 B = mesh.point(vB);
    Kernel::Point_3 C1 = mesh.point(vC1);
    Kernel::Point_3 C2 = mesh.point(vC2);
    Kernel::Vector_3 up(0.0, 1.0, 0.0);
    Kernel::Vector_3 x = B - A;
    Kernel::Vector_3 m = CGAL::cross_product(x, up);
    double s1 = CGAL::to_double(m * (C1 - A));
    double s2 = CGAL::to_double(m * (C2 - A));
    // class C
    if (s1 < 0.0 && s2 < 0.0)
        return 'C';

    // class D
    if (s1 > 0.0 && s2 > 0.0)
        return 'D';

    return 'X'; // something error
}


int compute_class_sum(
    const Surface_mesh& mesh,
    const std::set<Surface_mesh::Edge_index>& L,
    const std::set<Surface_mesh::Edge_index>& diagonal_edges)
{
    int l = 0;
    int cA = 0, cB = 0, cC = 0, cD = 0;

    for (auto e : mesh.edges()) 
    {
        // diagonal
        if (diagonal_edges.count(e))
            continue;

        // L
        if (L.count(e)) {
            l++;
            continue;
        }

        char cl = classify_edge_class(mesh, e);

        if      (cl == 'A') cA++;
        else if (cl == 'B') cB++;
        else if (cl == 'C') cC++;
        else if (cl == 'D') cD++;
        else {
            std::cout << "classified error\n";
            continue;
        } // X is error
    }

    int maxc = std::max({cA, cB, cC, cD});
    if (cA == maxc) cA = 0;
    else if (cB == maxc) cB = 0;
    else if (cC == maxc) cC = 0;
    else if (cD == maxc) cD = 0;

    return cA + cB + cC + cD + l;
}



double evaluate_gap(
    const Surface_mesh& mesh,
    const std::set<Surface_mesh::Edge_index>& L,
    const std::set<Surface_mesh::Edge_index>& diagonal_edges)
{
    int m = 0;

    // count edges without diagonal edges
    for (auto e : mesh.edges()) {
        if (!diagonal_edges.count(e))
            m++;
    }

    // check score of 4 class
    int score = compute_class_sum(mesh, L, diagonal_edges);

    // 5m/6
    double boundary = (5.0 * m) / 6.0;

    // gap
    double gap = boundary - score;

    std::cout << "[Result]\n";
    std::cout << "m = " << m << "\n";
    std::cout << "score(L) = " << score << "\n";
    std::cout << "5m/6 = " << boundary << "\n";
    std::cout << "gap = (5m/6 - score) = " << gap << "\n";
    if(gap > 0) std::cout << "true\n";
    else std::cout << "false\n";
    return gap;
}

void write_csv(const std::string& filename, int m, int score, double boundary, double gap)
{
    std::ofstream file(filename, std::ios::app);
    file << m << "," << score << "," << boundary << "," << gap << "\n";
}

struct AABB {
    double xmin, xmax, ymin, ymax, zmin, zmax;
};

bool aabb_inside(const AABB& inner, const AABB& outer) {
    return inner.xmin >= outer.xmin && inner.xmax <= outer.xmax &&
           inner.ymin >= outer.ymin && inner.ymax <= outer.ymax &&
           inner.zmin >= outer.zmin && inner.zmax <= outer.zmax;
}


int main() {
    
    int valid = 0; // counter of loop
    const int TARGET1 = 1; // Number of loop
    const int TARGET2 = 50; // Number of combinations
    int trials1 = 0; // counter of loop try
    const int MAX_TRIALS = 1000000000;

    Surface_mesh final_mesh;
    std::map<Surface_mesh::Face_index, Color> face_colors;
    std::map<Surface_mesh::Edge_index, Color> edge_colors;
    std::set<Surface_mesh::Edge_index> reflex_edges;
    std::set<Surface_mesh::Edge_index> diagonal_edges;
    std::set<Surface_mesh::Edge_index> guard_edges;

    bool have_mesh_to_draw = false;

    
    if (!glfwInit())
        return -1;

    GLFWwindow* window = glfwCreateWindow(1600, 1200, "CGAL Viewer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    const GLubyte* renderer = glGetString(GL_RENDERER);
    const GLubyte* version = glGetString(GL_VERSION);
    std::cout << "Renderer: " << renderer << "\n";
    std::cout << "OpenGL version supported: " << version << "\n";
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);


    // create csv file
    std::string csv_name = "Test_name.csv";
    std::ofstream file(csv_name);
    // header
    file << "m,score,5m/6,gap\n"; 

    std::srand(static_cast<unsigned>(time(NULL)));

    while(valid < TARGET1 && trials1 < MAX_TRIALS){
        trials1++;

        final_mesh.clear();
        face_colors.clear();
        edge_colors.clear();
        reflex_edges.clear();
        diagonal_edges.clear();
        guard_edges.clear();

        int success = 0;
        int trials2  = 0; // counter of conbine trying
        //int flag = 1; // flag for the polyhedron to use
        int flag = rand()%5; // if random
        // counter of polyhedron
        int tetra = 0;
        int cube = 0; 
        int octa = 0;
        int dodeca = 0;
        int icosa = 0;
        int er = 0;

        Nef_polyhedron polyhedron;
        Nef_polyhedron New_polyhedron;

        // create random polyhedron
        double c1 = rand()%20+1;
        double c2 = rand()%20+1;
        double c3 = rand()%20+1;
        double i1 = (rand()%2==0) ? -1.0 : 1.0;
        double i2 = (rand()%2==0) ? -1.0 : 1.0;
        double i3 = (rand()%2==0) ? -1.0 : 1.0;
        double size = rand()%30+5;

        Point first_center(
            Kernel::FT(c1 * i1),
            Kernel::FT(c2 * i2),
            Kernel::FT(c3 * i3)
        );

        if(flag == 0){ // conbine tetrahedra
            polyhedron = make_tetrahedron(Point(c1*i1, c2*i2, c3*i3), size);
            tetra++;
        }else if(flag == 1){ //conbine cube
            polyhedron = make_cube(Point(c1*i1, c2*i2, c3*i3), size);
            cube++;
        }else if(flag == 2){ // conbine octahedra
            polyhedron = make_Octahedron(Point(c1*i1, c2*i2, c3*i3), size);
            octa++;
        }
        else if(flag == 3){ // conbine dodecahedra
            polyhedron = make_Dodecahedron(Point(c1*i1, c2*i2, c3*i3), size);
            dodeca++;
        }else if(flag == 4){ // conbine icosahedra
            polyhedron = make_Icosahedron(Point(c1*i1, c2*i2, c3*i3), size);
            icosa++;
        }else{ // if error, conbine cube
            polyhedron = make_cube(Point(c1*i1, c2*i2, c3*i3), size);
            cube++;
            er++;
        }
        double r0 = size * 2.0;
        AABB current_box{
        CGAL::to_double(first_center.x()) - r0,
        CGAL::to_double(first_center.x()) + r0,
        CGAL::to_double(first_center.y()) - r0,
        CGAL::to_double(first_center.y()) + r0,
        CGAL::to_double(first_center.z()) - r0,
        CGAL::to_double(first_center.z()) + r0
        };

        success++;

        while (success < TARGET2 && trials2 < MAX_TRIALS) {
            trials2++;
            flag = rand()%5; // flag Update if random

            c1 = rand()%(20+(2*success))+1;
            c2 = rand()%(20+(2*success))+1;
            c3 = rand()%(20+(2*success))+1;
            i1 = (rand()%2==0) ? -1.0 : 1.0;
            i2 = (rand()%2==0) ? -1.0 : 1.0;
            i3 = (rand()%2==0) ? -1.0 : 1.0;
            size = rand()%30+5;

            if(flag == 0) New_polyhedron = make_tetrahedron(Point(c1*i1, c2*i2, c3*i3), size);
            else if(flag == 1) New_polyhedron = make_cube(Point(c1*i1, c2*i2, c3*i3), size);
            else if(flag == 2) New_polyhedron = make_Octahedron(Point(c1*i1, c2*i2, c3*i3), size);
            else if(flag == 3) New_polyhedron = make_Dodecahedron(Point(c1*i1, c2*i2, c3*i3), size);
            else if(flag == 4) New_polyhedron = make_Icosahedron(Point(c1*i1, c2*i2, c3*i3), size);
            else New_polyhedron = make_cube(Point(c1*i1, c2*i2, c3*i3), size);


            Point new_center(
            Kernel::FT(c1 * i1),
            Kernel::FT(c2 * i2),
            Kernel::FT(c3 * i3)
            );

            double r = size * 2.0;

            AABB new_box{
            CGAL::to_double(new_center.x()) - r,
            CGAL::to_double(new_center.x()) + r,
            CGAL::to_double(new_center.y()) - r,
            CGAL::to_double(new_center.y()) + r,
            CGAL::to_double(new_center.z()) - r,
            CGAL::to_double(new_center.z()) + r
            };
            

            if (aabb_inside(new_box, current_box)){
                continue;
            }


            if (!polyhedron.intersection(New_polyhedron).is_empty()) {

                try {
                    Nef_polyhedron tmp = polyhedron + New_polyhedron;

                    if (tmp.is_simple()){
                    polyhedron = tmp;
                    success++;
                    current_box.xmin = std::min(current_box.xmin, new_box.xmin);
                    current_box.xmax = std::max(current_box.xmax, new_box.xmax);
                    current_box.ymin = std::min(current_box.ymin, new_box.ymin);
                    current_box.ymax = std::max(current_box.ymax, new_box.ymax);
                    current_box.zmin = std::min(current_box.zmin, new_box.zmin);
                    current_box.zmax = std::max(current_box.zmax, new_box.zmax);
                    if(flag == 0) tetra++;
                    else if(flag == 1) cube++;
                    else if(flag == 2) octa++;
                    else if(flag == 3) dodeca++;
                    else if(flag == 4) icosa++;
                    else{ cube++;
                        er++;
                    }
                    continue;
                    }
                
                } catch (...) {

                    continue;
                }
            }
        }
        
        if(success < TARGET2){
            std::cout << "Warning: only " << success << " polyhedra merged\n";
            continue;
        }

        std::cout << "Tetrahedra: " << tetra << " merged\n";
        std::cout << "Cube: " << cube << " merged\n";
        std::cout << "Octahedra: " << octa << " merged\n";
        std::cout << "Icosahedra: " << icosa << " merged\n";
        std::cout << "Dodecahedra: " << dodeca << " merged\n";
        std::cout << "error:" << er << "\n"; 


        CGAL::convert_nef_polyhedron_to_polygon_mesh(polyhedron, final_mesh);

        // Added this 2 lines, for making triangulation and drawing certainly.
        assign_random_face_colors(final_mesh);
        CGAL::Polygon_mesh_processing::triangulate_faces(final_mesh);

        triangulate_faces_with_colors(final_mesh);


        if (!PMP::is_outward_oriented(final_mesh)) {
        PMP::reverse_face_orientations(final_mesh);
        }
       
        
        

        try {
        CGAL::convert_nef_polyhedron_to_polygon_mesh(polyhedron, final_mesh);
        } catch (...) {
            std::cout << "convert failed\n";
            continue;
        }

        if (!CGAL::is_triangle_mesh(final_mesh)) {
            try {
                PMP::triangulate_faces(final_mesh);
                triangulate_faces_with_colors(final_mesh);
            } catch (...) {
                std::cout << "triangulate failed\n";
                continue;
            }
        }

        if (!PMP::is_outward_oriented(final_mesh)) {
            PMP::reverse_face_orientations(final_mesh);
        }


        classify_edges_ver3(final_mesh, edge_colors, reflex_edges, diagonal_edges);

        auto graph = One_skelton(final_mesh);
        guard_edges = GreedyEdgeGuards(final_mesh, diagonal_edges);

        double gap = evaluate_gap(final_mesh, guard_edges, diagonal_edges);
    
        int m = 0;
        for (auto e : final_mesh.edges())
            if (!diagonal_edges.count(e)) m++;

        int score = compute_class_sum(final_mesh, guard_edges, diagonal_edges);
        double boundary = (5.0 * m) / 6.0;

        if (m < 1000){
            continue;
        }

        valid++;

        // save to CSV
        write_csv(csv_name, m, score, boundary, gap);

        for (auto face : final_mesh.faces()) {
            face_colors[face] = {
                static_cast<float>(rand()) / RAND_MAX,
                static_cast<float>(rand()) / RAND_MAX,
                static_cast<float>(rand()) / RAND_MAX
            };
        }
        have_mesh_to_draw = true;
    }
    

    if(have_mesh_to_draw){

        while (!glfwWindowShouldClose(window)) {
            int width, height;
            glfwGetFramebufferSize(window, &width, &height);
            float aspect = (float)width / (float)height;

            glViewport(0, 0, width, height);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            set_perspective(45.0f, aspect, 0.1f, 10000.0f);

            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glTranslatef(0.0f, 0.0f, -1500.0f);
            glRotatef(rotation_x, 1.0f, 0.0f, 0.0f);
            glRotatef(rotation_y, 0.0f, 1.0f, 0.0f);

            // left:final_mesh
            glPushMatrix();
            glTranslatef(-300, 0.0f, 0.0f);  // align to the left
            draw_mesh(final_mesh, face_colors, edge_colors, reflex_edges);
            glPopMatrix();

            // right:skelton + guard
            glPushMatrix();
            glTranslatef(+300, 0.0f, 0.0f);  // align to the right
            draw_skeleton(final_mesh);    // 1-skelton
            draw_guard_edges(final_mesh, guard_edges);    // edge guards(yellow)
            glPopMatrix();


            glfwSwapBuffers(window);
            glfwPollEvents();
        }
    }
    glfwDestroyWindow(window);
    glfwTerminate();    

    return 0;
}
