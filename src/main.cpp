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

Nef_polyhedron make_cube(Point center, double size) {
    double s = size / 2.0;
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

Nef_polyhedron make_tetrahedron(Point center, double size) {
    double s = size / 2.0;
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

Nef_polyhedron make_random_tetrahedron(Point center, double size) {
    int s = size;
    while(true){
        std::srand(time(NULL));
        double s1_1 = rand()%s;
        double s1_2 = rand()%s;
        double s1_3 = rand()%s;
        double i1_1 = (rand()%2==0) ? -1.0 : 1.0;
        double i1_2 = (rand()%2==0) ? -1.0 : 1.0;
        double i1_3 = (rand()%2==0) ? -1.0 : 1.0;
        double s2_1 = rand()%s;
        double s2_2 = rand()%s;
        double s2_3 = rand()%s;
        double i2_1 = (rand()%2==0) ? -1.0 : 1.0;
        double i2_2 = (rand()%2==0) ? -1.0 : 1.0;
        double i2_3 = (rand()%2==0) ? -1.0 : 1.0;
        double s3_1 = rand()%s;
        double s3_2 = rand()%s;
        double s3_3 = rand()%s;
        double i3_1 = (rand()%2==0) ? -1.0 : 1.0;
        double i3_2 = (rand()%2==0) ? -1.0 : 1.0;
        double i3_3 = (rand()%2==0) ? -1.0 : 1.0;
        double s4_1 = rand()%s;
        double s4_2 = rand()%s;
        double s4_3 = rand()%s;
        double i4_1 = (rand()%2==0) ? -1.0 : 1.0;
        double i4_2 = (rand()%2==0) ? -1.0 : 1.0;
        double i4_3 = (rand()%2==0) ? -1.0 : 1.0;
        std::vector<Point> points = {
            Point(center.x()+s1_1*i1_1, center.y()+s1_2*i1_2, center.z()+s1_3*i1_3),
            Point(center.x()+s2_1*i2_1, center.y()+s2_2*i2_2, center.z()+s2_3*i2_3),
            Point(center.x()+s3_1*i3_1, center.y()+s3_2*i3_2, center.z()+s3_3*i3_3),
            Point(center.x()+s4_1*i4_1, center.y()+s4_2*i4_2, center.z()+s4_3*i4_3),
        };

        if (is_degenerate_tetrahedron(points[0], points[1], points[2], points[3])) {
            continue; // 再生成
        }

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
}

Nef_polyhedron make_random_tetrahedron2(Point center, double size) {
    while (true) {
        std::vector<Point> points;
        for (int i = 0; i < 4; ++i) {
            double x = ((rand() % 1000) / 1000.0 - 0.5) * size;
            double y = ((rand() % 1000) / 1000.0 - 0.5) * size;
            double z = ((rand() % 1000) / 1000.0 - 0.5) * size;
            points.emplace_back(center.x() + x, center.y() + y, center.z() + z);
        }

        if (is_degenerate_tetrahedron(points[0], points[1], points[2], points[3])) {
            continue; // 再生成
        }

        Surface_mesh mesh;
        std::vector<Surface_mesh::Vertex_index> vi;
        for (int i = 0; i < 4; ++i) {
            vi.push_back(mesh.add_vertex(points[i]));
        }

        mesh.add_face(vi[0], vi[2], vi[1]);
        mesh.add_face(vi[0], vi[1], vi[3]);
        mesh.add_face(vi[0], vi[3], vi[2]);
        mesh.add_face(vi[1], vi[2], vi[3]);

        return Nef_polyhedron(mesh);
    }
}
Nef_polyhedron make_Dodecahedron(Point center, double size){
    double s = size / 2.0;
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
    double s = size / 2.0;
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

Nef_polyhedron make_Octahedron(Point center, double size) {
    double s = size;

    // 6つの頂点（Z軸に上下、XY平面に正方形の4点）
    std::vector<Point> points = {
        Point(center.x(),     center.y(),     center.z() + s), // 上 (頂点0)
        Point(center.x(),     center.y(),     center.z() - s), // 下 (頂点1)
        Point(center.x() + s, center.y(),     center.z()),     // 右 (頂点2)
        Point(center.x() - s, center.y(),     center.z()),     // 左 (頂点3)
        Point(center.x(),     center.y() + s, center.z()),     // 前 (頂点4)
        Point(center.x(),     center.y() - s, center.z())      // 後 (頂点5)
    };

    Surface_mesh mesh;
    auto v = [&](int i) { return mesh.add_vertex(points[i]); };
    std::vector<Surface_mesh::Vertex_index> vi(6);
    for (int i = 0; i < 6; ++i) vi[i] = v(i);

    auto add_face = [&](int a, int b, int c) {
        mesh.add_face(vi[a], vi[b], vi[c]);
    };

    // 上半球 4面
    add_face(0, 2, 4);
    add_face(0, 4, 3);
    add_face(0, 3, 5);
    add_face(0, 5, 2);

    // 下半球 4面
    add_face(1, 4, 2);
    add_face(1, 3, 4);
    add_face(1, 5, 3);
    add_face(1, 2, 5);

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
    const std::set<Surface_mesh::Edge_index>& reflex_edges_set) {
    // draw faces
    for (auto face : mesh.faces()) {
        glBegin(GL_TRIANGLE_FAN);
        auto color = face_colors.at(face);
        glColor3f(color.r, color.g, color.b);

        for (auto v : CGAL::vertices_around_face(mesh.halfedge(face), mesh)) {
            Point p = mesh.point(v);
            glVertex3f(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
        }
        glEnd();
    }

    // draw edges
    glLineWidth(3.0f);
    glBegin(GL_LINES);
    for (auto edge : mesh.edges()) {
        auto color = edge_colors.at(edge);   // ← クラスA〜Dごとの色を取得
        glColor3f(color.r, color.g, color.b);

        // reflex 辺っぽい色なら太く描く
        // （classify_edges_ver3 では淡色＝reflex なので、g成分などで判別可）
        bool isReflexEdge = reflex_edges_set.find(edge) != reflex_edges_set.end();

    
        auto h = mesh.halfedge(edge);
        Point p1 = mesh.point(mesh.source(h));
        Point p2 = mesh.point(mesh.target(h));

        // Reflex辺は見た目だけ手前に少しずらす
        if (isReflexEdge) {
            auto f1 = mesh.face(h);
            auto f2 = mesh.face(mesh.opposite(h));

            Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
            Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);
            Kernel::Vector_3 avg_n = (n1 + n2) / 2.0;

            double len = std::sqrt(CGAL::to_double(avg_n.squared_length()));
            if (len > 1e-8) avg_n = avg_n / len;

            double offset = 0.1;  // 👈 この値で“浮かせる量”を調整（大きすぎるとズレる）
            p1 = Point(p1.x() + avg_n.x() * offset,
                       p1.y() + avg_n.y() * offset,
                       p1.z() + avg_n.z() * offset);
            p2 = Point(p2.x() + avg_n.x() * offset,
                       p2.y() + avg_n.y() * offset,
                       p2.z() + avg_n.z() * offset);
        }

        glLineWidth(isReflexEdge ? 8.0f : 3.0f); // reflexは太く
        glVertex3f(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p1.z()));
        glVertex3f(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), CGAL::to_double(p2.z()));
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
void One_skelton(){

}

// greedy algorithm
void GreedyAlgorithm(){

}

bool is_same_polygon_face(const Surface_mesh& mesh, Surface_mesh::Face_index f1, Surface_mesh::Face_index f2) {
    if (f1 == Surface_mesh::null_face() || f2 == Surface_mesh::null_face())
        return false;

    Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
    Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

    // 法線方向がほぼ同じ（1.0に近い）
    return std::fabs(CGAL::to_double(n1 * n2 - 1.0)) < 1e-6;
}

// エッジ分類関数
void classify_edges(const Surface_mesh& mesh) {
    int real_edges = 0;
    int reflex_edges = 0;
    int convex_edges = 0;
    int diagonal_edges = 0;

    for (auto edge : mesh.edges()) {
        auto h1 = mesh.halfedge(edge, 0);

        // 境界エッジは実エッジとしてカウント
        if (mesh.is_border(h1)) {
            real_edges++;
            continue;
        }

        auto h2 = mesh.opposite(h1);
        auto f1 = mesh.face(h1);
        auto f2 = mesh.face(h2);

        // 同一平面上の面なら対角線
        if (is_same_polygon_face(mesh, f1, f2)) {
            diagonal_edges++;
            continue;
        }

        // 面の法線ベクトル
        Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
        Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

        // 判定用のベクトル（h1の始点からh2の次の頂点へ）
        auto v_source = mesh.source(h1);
        auto v_other = mesh.target(mesh.next(h2));
        Kernel::Vector_3 test_vec = mesh.point(v_other) - mesh.point(v_source);

        // 凹凸判定（右ねじ方向で外側なら凸、内側なら凹）
        if (test_vec * n1 > 0) {
            reflex_edges++;
        } else {
            convex_edges++;
        }

        real_edges++;
    }

    std::cout << "[Edge Classification Result]\n";
    std::cout << "Total Edges    : " << mesh.number_of_edges() << "\n";
    std::cout << "Real Edges     : " << real_edges << "\n";
    std::cout << "Reflex Edges   : " << reflex_edges << "\n";
    std::cout << "Convex Edges   : " << convex_edges << "\n";
    std::cout << "Diagonal Edges : " << diagonal_edges << "\n";
}


void classify_edges_ver2(const Surface_mesh& mesh,
                         std::map<Surface_mesh::Edge_index, Color>& edge_colors) {
    int real_edges = 0;
    int reflex_edges = 0;
    int convex_edges = 0;
    int diagonal_edges = 0;

    int classA = 0, classB = 0, classC = 0, classD = 0;

    for (auto edge : mesh.edges()) {
        auto h1 = mesh.halfedge(edge, 0);

        // 境界エッジは実エッジとしてカウントしてスキップ
        if (mesh.is_border(h1)) {
            real_edges++;
            edge_colors[edge] = {1.0f, 1.0f, 1.0f}; // 白
            continue;
        }

        auto h2 = mesh.opposite(h1);
        auto f1 = mesh.face(h1);
        auto f2 = mesh.face(h2);

        // --- 面の法線 ---
        Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
        Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

        // --- 法線の内積による平面判定 ---
        double dot = CGAL::to_double(n1 * n2);
        if (std::fabs(dot - 1.0) < 1e-6) {
            // 同一平面上 → 対角線として扱う
            diagonal_edges++;
            edge_colors[edge] = {0.7f, 0.7f, 0.7f}; // 灰色
            continue;
        }

        // --- 凹凸判定 ---
        auto v_source = mesh.source(h1);
        auto v_other = mesh.target(mesh.next(h2));
        Kernel::Vector_3 test_vec = mesh.point(v_other) - mesh.point(v_source);

        if (test_vec * n1 > 0)
            reflex_edges++;
        else
            convex_edges++;

        real_edges++;

        // --- クラス分類 ---
        double x1 = CGAL::to_double(n1.x());
        double y1 = CGAL::to_double(n1.y());
        double x2 = CGAL::to_double(n2.x());
        double y2 = CGAL::to_double(n2.y());

        bool isA = (y1 >= 0 && y2 >= 0);
        bool isB = (y1 < 0 && y2 < 0);
        bool isC = (
            (x1 >= 0 && x2 >= 0) ||
            ((x1 >= 0 && y1 >= 0) && (x2 < 0 && y2 < 0)) ||
            ((x2 >= 0 && y2 >= 0) && (x1 < 0 && y1 < 0))
        );
        bool isD = (
            (x1 < 0 && x2 < 0) ||
            ((x1 < 0 && y1 >= 0) && (x2 >= 0 && y2 < 0)) ||
            ((x2 < 0 && y2 >= 0) && (x1 >= 0 && y1 < 0))
        );

        if (isA) { classA++; edge_colors[edge] = {1.0f, 0.0f, 0.0f}; } // 赤
        else if (isB) { classB++; edge_colors[edge] = {0.0f, 0.0f, 1.0f}; } // 青
        else if (isC) { classC++; edge_colors[edge] = {0.0f, 1.0f, 0.0f}; } // 緑
        else if (isD) { classD++; edge_colors[edge] = {1.0f, 1.0f, 0.0f}; } // 黄
        else { edge_colors[edge] = {1.0f, 1.0f, 1.0f}; } // デフォルト白
    }

    std::cout << "[Edge Classification Result]\n";
    std::cout << "Total Edges    : " << mesh.number_of_edges() << "\n";
    std::cout << "Real Edges     : " << real_edges << "\n";
    std::cout << "Reflex Edges   : " << reflex_edges << "\n";
    std::cout << "Convex Edges   : " << convex_edges << "\n";
    std::cout << "Diagonal Edges : " << diagonal_edges << "\n\n";

    std::cout << "[A–D Class Result]\n";
    std::cout << "Class A (Y>=0 both): " << classA << "\n";
    std::cout << "Class B (Y<0 both):  " << classB << "\n";
    std::cout << "Class C (X>=0 or mixed +Y/-Y): " << classC << "\n";
    std::cout << "Class D (X<0 or mixed +Y/-Y):  " << classD << "\n";
}

void classify_edges_ver3(const Surface_mesh& mesh,
                         std::map<Surface_mesh::Edge_index, Color>& edge_colors,
                         std::set<Surface_mesh::Edge_index>& reflex_edges_set) {
    int real_edges = 0;
    int reflex_edges = 0;
    int convex_edges = 0;
    int diagonal_edges = 0;

    // 各クラスのカウント
    int classA_convex = 0, classA_reflex = 0;
    int classB_convex = 0, classB_reflex = 0;
    int classC_convex = 0, classC_reflex = 0;
    int classD_convex = 0, classD_reflex = 0;

    for (auto edge : mesh.edges()) {
        auto h1 = mesh.halfedge(edge, 0);

        // 境界エッジはスキップ
        if (mesh.is_border(h1)) {
            real_edges++;
            edge_colors[edge] = {1.0f, 1.0f, 1.0f}; // 白
            continue;
        }

        auto h2 = mesh.opposite(h1);
        auto f1 = mesh.face(h1);
        auto f2 = mesh.face(h2);

        // 面の法線ベクトル
        Kernel::Vector_3 n1 = PMP::compute_face_normal(f1, mesh);
        Kernel::Vector_3 n2 = PMP::compute_face_normal(f2, mesh);

        // 平面が同じなら対角線
        double dot = CGAL::to_double(n1 * n2);
        if (std::fabs(dot - 1.0) < 1e-6) {
            diagonal_edges++;
            edge_colors[edge] = {0.7f, 0.7f, 0.7f};
            continue;
        }

        // --- 凹凸判定 ---
        auto v_source = mesh.source(h1);
        auto v_other = mesh.target(mesh.next(h2));
        Kernel::Vector_3 test_vec = mesh.point(v_other) - mesh.point(v_source);

        bool isReflex = (test_vec * n1 > 0);
        if (isReflex){
            reflex_edges++;
            reflex_edges_set.insert(edge);  // ← reflexエッジを登録
        }
        else
            convex_edges++;
        real_edges++;

        // --- A～D クラス分類条件 ---
        double x1 = CGAL::to_double(n1.x());
        double y1 = CGAL::to_double(n1.y());
        double x2 = CGAL::to_double(n2.x());
        double y2 = CGAL::to_double(n2.y());

        bool isA = (y1 >= 0 && y2 >= 0);
        bool isB = (y1 < 0 && y2 < 0);
        bool isC = (
            (x1 >= 0 && x2 >= 0) ||
            ((x1 >= 0 && y1 >= 0) && (x2 < 0 && y2 < 0)) ||
            ((x2 >= 0 && y2 >= 0) && (x1 < 0 && y1 < 0))
        );
        bool isD = (
            (x1 < 0 && x2 < 0) ||
            ((x1 < 0 && y1 >= 0) && (x2 >= 0 && y2 < 0)) ||
            ((x2 < 0 && y2 >= 0) && (x1 >= 0 && y1 < 0))
        );

        // --- 色分けとカウント ---
        if (isA) {
            if (isReflex) { classA_reflex++; edge_colors[edge] = {1.0f, 0.5f, 0.5f}; } // Reflex A2（淡赤）
            else          { classA_convex++; edge_colors[edge] = {1.0f, 0.0f, 0.0f}; } // Convex A1（赤）
        }
        else if (isB) {
            if (isReflex) { classB_reflex++; edge_colors[edge] = {0.5f, 0.5f, 1.0f}; } // Reflex B2（淡青）
            else          { classB_convex++; edge_colors[edge] = {0.0f, 0.0f, 1.0f}; } // Convex B1（青）
        }
        else if (isC) {
            if (isReflex) { classC_reflex++; edge_colors[edge] = {0.5f, 1.0f, 0.5f}; } // Reflex C2（淡緑）
            else          { classC_convex++; edge_colors[edge] = {0.0f, 1.0f, 0.0f}; } // Convex C1（緑）
        }
        else if (isD) {
            if (isReflex) { classD_reflex++; edge_colors[edge] = {1.0f, 1.0f, 0.5f}; } // Reflex D2（淡黄）
            else          { classD_convex++; edge_colors[edge] = {1.0f, 1.0f, 0.0f}; } // Convex D1（黄）
        }
        else {
            edge_colors[edge] = {1.0f, 1.0f, 1.0f};
        }
    }

    // --- 結果表示 ---
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
    std::cout << "D2 (X<0 or mixed +Y/-Y) Reflex : " << classD_reflex << "\n";
}





int main() {
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

    // create two cubes and merge them
    //Nef_polyhedron cube1 = make_cube(Point(0.25, 0.25, 0.25), 1.0);
    //Nef_polyhedron cube2 = make_cube(Point(-0.25, -0.25, -0.25), 1.0);
    //Nef_polyhedron cubes = cube1 + cube2;

    // create regular tetrahydron
    //Nef_polyhedron tetrahedra = make_tetrahedron(Point(0.25, 0.25, 0.25), 1.0);

    // create random tetrahydron
    //Nef_polyhedron tetrahedra = make_random_tetrahedron(Point(0.25, 0.25, 0.25), 2.0);
    //for(int i=0 ; i<8 ; i++){
    //    Nef_polyhedron New_tetrahedra = make_random_tetrahedron(Point(0.25, 0.25, 0.25), 3.0*(i+2));
    //    if(!tetrahedra.intersection(New_tetrahedra).is_empty()){
    //        tetrahedra = tetrahedra + New_tetrahedra;
    //    }
    //}

    // create regular tetrahydron

    //Nef_polyhedron tetrahedra = make_Icosahedron(Point(0.25, 0.25, 0.25), 20.0);
    //Nef_polyhedron tetrahedra = make_Dodecahedron(Point(0.25, 0.25, 0.25), 20.0);
    //Nef_polyhedron tetrahedra = make_Octahedron(Point(0.25, 0.25, 0.25), 20.0);

    // create random tetrahydron
    std::srand(time(NULL));
    double c1 = rand()%20+1;
    double c2 = rand()%20+1;
    double c3 = rand()%20+1;
    double i1 = (rand()%2==0) ? -1.0 : 1.0;
    double i2 = (rand()%2==0) ? -1.0 : 1.0;
    double i3 = (rand()%2==0) ? -1.0 : 1.0;
    double size = rand()%30+1;
    Nef_polyhedron tetrahedra = make_Icosahedron(Point(c1*i1, c2*i2, c3*i3), size);
    //Nef_polyhedron tetrahedra = make_Dodecahedron(Point(c1*i1, c2*i2, c3*i3), size);
    //Nef_polyhedron tetrahedra = make_Octahedron(Point(c1*i1, c2*i2, c3*i3), size);
    for(int i=0 ; i<100 ; i++){
        c1 = rand()%(20+i)+1;
        c2 = rand()%(20+i)+1;
        c3 = rand()%(20+i)+1;
        i1 = (rand()%2==0) ? -1.0 : 1.0;
        i2 = (rand()%2==0) ? -1.0 : 1.0;
        i3 = (rand()%2==0) ? -1.0 : 1.0;
        size = rand()%30+1;
        Nef_polyhedron New_tetrahedra = make_Icosahedron(Point(c1*i1, c2*i2, c3*i3), size);
        //Nef_polyhedron New_tetrahedra = make_Dodecahedron(Point(c1*i1, c2*i2, c3*i3), size);
        //Nef_polyhedron New_tetrahedra = make_Octahedron(Point(c1*i1, c2*i2, c3*i3), size);
        if(!tetrahedra.intersection(New_tetrahedra).is_empty()){
            tetrahedra = tetrahedra + New_tetrahedra;
        }
    }

    //Nef_polyhedron cube3 = make_cube(Point(1000.75, 0, 0), 2000.0);
    //Kernel::Vector_3 axis(0, 0, 1); // rotate around z-axis

    //Nef_polyhedron rotated = rotate_nef_polyhedron(cube3, axis, 45);
    //Nef_polyhedron rotated = rotate_nef_polyhedron(tetrahedra1, axis, 45);

    //Nef_polyhedron result = cubes.intersection(rotated.complement());
    //Nef_polyhedron result = tetrahedra1.intersection(rotated.complement());

    //rotated = rotate_nef_polyhedron(rotated, axis, 180);
    //result = result.intersection(rotated.complement());

    Surface_mesh final_mesh;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(tetrahedra, final_mesh);

    // Added this 2 lines, for making triangulation and drawing certainly.
    assign_random_face_colors(final_mesh);
    CGAL::Polygon_mesh_processing::triangulate_faces(final_mesh);

    triangulate_faces_with_colors(final_mesh);


    if (!PMP::is_outward_oriented(final_mesh)) {
    PMP::reverse_face_orientations(final_mesh);
    }

    std::map<Surface_mesh::Edge_index, Color> edge_colors;
    std::set<Surface_mesh::Edge_index> reflex_edges;

    classify_edges_ver3(final_mesh, edge_colors, reflex_edges);


    std::map<Surface_mesh::Face_index, Color> face_colors;
    for (auto face : final_mesh.faces()) {
        face_colors[face] = {
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX,
            static_cast<float>(rand()) / RAND_MAX
        };
    }



    while (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float aspect = (float)width / (float)height;

        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        set_perspective(45.0f, aspect, 0.1f, 1000.0f); //描画距離

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -200.0f); //描画視点をどれぐらい引くかのパラメータ
        glRotatef(rotation_x, 1.0f, 0.0f, 0.0f);
        glRotatef(rotation_y, 0.0f, 1.0f, 0.0f);

        draw_mesh(final_mesh, face_colors, edge_colors, reflex_edges);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
