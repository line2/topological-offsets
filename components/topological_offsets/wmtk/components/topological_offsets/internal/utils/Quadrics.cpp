#include "Quadrics.hpp"

#include <cstdint>
#include <wmtk/simplex/faces_single_dimension.hpp>
#include <wmtk/simplex/top_dimension_cofaces.hpp>
#include <wmtk/simplex/top_dimension_cofaces_iterable.hpp>

namespace wmtk::components::internal::utils {

namespace {

struct Sample
{
    SimpleBVH::VectorMax3d point;
    SimpleBVH::VectorMax3d normal;
    double offset_distance;
    double distance;
    double weight;
    double area;
};

Sample get_point_sample(const OffsetBvh& bvh, const SimpleBVH::VectorMax3d& p)
{
    const auto [sq_dist, offset_distance, nearest_point] =
        bvh.sq_dist_offset_distance_nearest_point(p);

    Sample sample;
    sample.point = nearest_point;
    sample.normal = (p - nearest_point).normalized(); // TODO use distance to normalize
    sample.offset_distance = offset_distance;
    sample.distance = (p - nearest_point).norm();
    return sample;
}

void get_triangle_samples_and_area(
    const OffsetBvh& bvh,
    const attribute::MeshAttributeHandle& pos_handle,
    const simplex::Simplex& s,
    std::vector<Sample>& samples)
{
    const Mesh& m = pos_handle.mesh();
    const auto pos_acc = m.create_const_accessor(pos_handle.as<double>());

    const auto vertices = simplex::faces_single_dimension_tuples(m, s, PrimitiveType::Vertex);
    assert(vertices.size() == 3);

    const Eigen::Vector3d p0 = pos_acc.const_vector_attribute(vertices[0]);
    const Eigen::Vector3d p1 = pos_acc.const_vector_attribute(vertices[1]);
    const Eigen::Vector3d p2 = pos_acc.const_vector_attribute(vertices[2]);

    const Eigen::Vector3d p_mid = (p0 + p1 + p2) / 3.;

    const double area = (p1 - p0).cross(p2 - p0).norm();

    constexpr double u = 0.1;
    const std::array<Eigen::Vector3d, 4> sample_points{
        {p_mid, (1 - u) * p0 + u * p_mid, (1 - u) * p1 + u * p_mid, (1 - u) * p2 + u * p_mid}};

    std::vector<Sample> new_samples;
    new_samples.reserve(4);

    for (const Eigen::Vector3d& point : sample_points) {
        new_samples.emplace_back(get_point_sample(bvh, point));
    }

    new_samples[0].weight = 1;
    new_samples[1].weight = u;
    new_samples[2].weight = u;
    new_samples[3].weight = u;

    for (Sample& sam : new_samples) {
        sam.area = area;
    }

    std::copy(new_samples.begin(), new_samples.end(), std::back_inserter(samples));
}

void get_edge_samples_and_area(
    const OffsetBvh& bvh,
    const attribute::MeshAttributeHandle& pos_handle,
    const simplex::Simplex& s,
    std::vector<Sample>& samples)
{
    const Mesh& m = pos_handle.mesh();
    const auto pos_acc = m.create_const_accessor(pos_handle.as<double>());

    const auto vertices = simplex::faces_single_dimension_tuples(m, s, PrimitiveType::Vertex);
    assert(vertices.size() == 2);

    const Eigen::Vector2d p0 = pos_acc.const_vector_attribute(vertices[0]);
    const Eigen::Vector2d p1 = pos_acc.const_vector_attribute(vertices[1]);

    const Eigen::Vector2d p_mid = (p0 + p1) / 2.;

    const double area = (p1 - p0).norm();

    constexpr double u = 0.1;
    const std::array<Eigen::Vector2d, 3> sample_points{
        {p_mid, (1 - u) * p0 + u * p_mid, (1 - u) * p1 + u * p_mid}};

    std::vector<Sample> new_samples;
    new_samples.reserve(3);

    for (const Eigen::Vector2d& point : sample_points) {
        new_samples.emplace_back(get_point_sample(bvh, point));
    }

    new_samples[0].weight = 1;
    new_samples[1].weight = u;
    new_samples[2].weight = u;

    for (Sample& sam : new_samples) {
        sam.area = area;
    }

    std::copy(new_samples.begin(), new_samples.end(), std::back_inserter(samples));
}

} // namespace

Quadrics::Quadrics(const double a, const double b, const double c, const double d)
{
    Eigen::Vector4d nh(a, b, c, d);
    m_matrix = nh * nh.transpose();
}

Quadrics::Quadrics(const Eigen::Vector3d& p0, const Eigen::Vector3d& n)
{
    const double d = -n.transpose() * p0;
    Eigen::Vector4d nh(n[0], n[1], n[2], d);
    m_matrix = nh * nh.transpose();
}

Quadrics::Quadrics(const attribute::MeshAttributeHandle& quadric_handle, const simplex::Simplex& s)
{
    const Mesh& m = quadric_handle.mesh();
    assert(
        m.top_simplex_type() == PrimitiveType::Triangle ||
        m.top_simplex_type() == PrimitiveType::Edge);

    if (s.primitive_type() == PrimitiveType::Triangle) {
        m_matrix = Quadrics::attribute_to_matrix(quadric_handle, s.tuple());
    } else if (s.primitive_type() == PrimitiveType::Vertex) {
        m_matrix.setZero();
        auto cofaces = simplex::top_dimension_cofaces_iterable(m, s);
        for (const Tuple& f : cofaces) {
            const Eigen::Matrix4d m = Quadrics::attribute_to_matrix(quadric_handle, f);
            m_matrix += m;
        }
    } else {
        log_and_throw_error("Unknown primitive type for Quadrics.");
    }
}

Quadrics::Quadrics(
    const OffsetBvh& bvh,
    const attribute::MeshAttributeHandle& pos_handle,
    const attribute::MeshAttributeHandle& tet_tag_handle,
    const attribute::MeshAttributeHandle& bvh_mesh_position_handle,
    const simplex::Simplex& s)
{
    const Mesh& m = pos_handle.mesh();

    assert(!m.is_multi_mesh_root());
    assert(!m.has_child_mesh());
    assert(
        s.primitive_type() == PrimitiveType::Vertex || s.primitive_type() == PrimitiveType::Edge ||
        s.primitive_type() == PrimitiveType::Triangle); // only implemented for these

    m_matrix.setZero();

    if (s.primitive_type() == PrimitiveType::Triangle) {
        assert(m.top_simplex_type() == PrimitiveType::Triangle);

        std::vector<Sample> samples;
        get_triangle_samples_and_area(bvh, pos_handle, s, samples);

        Quadrics q(
            samples[0].point + samples[0].offset_distance * samples[0].normal,
            samples[0].normal);

        for (int64_t i = 1; i < samples.size(); ++i) {
            q += Quadrics(
                     samples[i].point + samples[i].offset_distance * samples[i].normal,
                     samples[i].normal) *
                 samples[i].weight;
        }

        q *= samples[0].area;

        m_matrix = q.m_matrix;

    } else if (s.primitive_type() == PrimitiveType::Edge) {
        assert(m.top_simplex_type() == PrimitiveType::Edge);
        assert(pos_handle.dimension() == 2);

        std::vector<Sample> samples;
        get_edge_samples_and_area(bvh, pos_handle, s, samples);

        Eigen::Vector3d q_point_3(samples[0].point[0], samples[0].point[1], 0);
        Eigen::Vector3d q_normal_3(samples[0].normal[0], samples[0].normal[1], 0);
        Quadrics q(q_point_3 + samples[0].offset_distance * q_normal_3, q_normal_3);

        for (int64_t i = 1; i < samples.size(); ++i) {
            Eigen::Vector3d q_point_3(samples[i].point[0], samples[i].point[1], 0);
            Eigen::Vector3d q_normal_3(samples[i].normal[0], samples[i].normal[1], 0);
            q += Quadrics(q_point_3 + samples[i].offset_distance * q_normal_3, q_normal_3) *
                 samples[i].weight;
        }

        q *= samples[0].area;

        m_matrix = q.m_matrix;
    } else if (s.primitive_type() == PrimitiveType::Vertex) {
        auto cofaces = simplex::top_dimension_cofaces_iterable(m, s);

        if (m.top_simplex_type() == PrimitiveType::Triangle) {
            assert(pos_handle.dimension() == 3);
            std::vector<Sample> samples;
            samples.reserve(6 * 4);

            for (const Tuple& f : cofaces) {
                const simplex::Simplex f_simplex(m, m.top_simplex_type(), f);
                get_triangle_samples_and_area(bvh, pos_handle, f_simplex, samples);
            }

            double dist_avg = 0;
            double area_sum = 0;
            double dist_min = std::numeric_limits<double>::max();
            for (const Sample& sam : samples) {
                dist_avg += sam.offset_distance * sam.area;
                area_sum += sam.area;
                dist_min = std::min(dist_min, sam.distance);
            }
            dist_avg /= area_sum;

            dist_avg = 0.5 * (dist_avg + dist_min);

            for (const Sample& sam : samples) {
                (*this) +=
                    Quadrics(sam.point + dist_avg * sam.normal, sam.normal) * sam.weight * sam.area;
            }
        } else if (m.top_simplex_type() == PrimitiveType::Edge) {
            assert(pos_handle.dimension() == 2);
            std::vector<Sample> samples;
            samples.reserve(2 * 3);

            for (const Tuple& f : cofaces) {
                const simplex::Simplex f_simplex(m, m.top_simplex_type(), f);
                get_edge_samples_and_area(bvh, pos_handle, f_simplex, samples);
            }

            double dist_avg = 0;
            double area_sum = 0;
            double dist_min = std::numeric_limits<double>::max();
            for (const Sample& sam : samples) {
                dist_avg += sam.offset_distance * sam.area;
                area_sum += sam.area;
                dist_min = std::min(dist_min, sam.distance);
            }
            dist_avg /= area_sum;

            dist_avg = 0.5 * (dist_avg + dist_min);

            for (const Sample& sam : samples) {
                Eigen::Vector3d q_point_3(sam.point[0], sam.point[1], 0);
                Eigen::Vector3d q_normal_3(sam.normal[0], sam.normal[1], 0);
                (*this) +=
                    Quadrics(q_point_3 + dist_avg * q_normal_3, q_normal_3) * sam.weight * sam.area;
            }
        } else {
            log_and_throw_error("Unknown mesh type in Quadrics.");
        }

        // auto cofaces = simplex::top_dimension_cofaces_iterable(m, s);
        // for (const Tuple& f : cofaces) {
        //     const simplex::Simplex f_simplex(m, m.top_simplex_type(), f);
        //     (*this) += Quadrics(
        //         bvh,
        //         pos_handle,
        //         tet_tag_handle,
        //         bvh_mesh_position_handle,
        //         f_simplex,
        //         offset_distance_handle,
        //         bvh_to_input_mesh);
        // }
    } else {
        log_and_throw_error("Unknown primitive type for Quadrics.");
    }
}

Quadrics::Quadrics(
    const SimpleBVH::BVH& bvh,
    const attribute::MeshAttributeHandle& pos_handle,
    const attribute::MeshAttributeHandle& tet_tag_handle,
    const attribute::MeshAttributeHandle& bvh_mesh_position_handle,
    const simplex::Simplex& s,
    const double offset_distance)
{
    assert(
        s.primitive_type() == PrimitiveType::Vertex || s.primitive_type() == PrimitiveType::Edge ||
        s.primitive_type() == PrimitiveType::Triangle); // only implemented for these

    const Mesh& m = pos_handle.mesh();

    m_matrix.setZero();

    if (s.primitive_type() == PrimitiveType::Triangle) {
        assert(!m.is_multi_mesh_root());
        assert(!m.has_child_mesh());
        assert(m.top_simplex_type() == PrimitiveType::Triangle);

        // const Mesh& embedding_mesh = tet_tag_handle.mesh();
        // const attribute::MeshAttributeHandle embedding_pos_handle =
        //     embedding_mesh.get_attribute_handle<double>("vertices", PrimitiveType::Vertex);

        const auto vertices = simplex::faces_single_dimension_tuples(m, s, PrimitiveType::Vertex);
        assert(vertices.size() == 3);

        const auto pos_acc = m.create_const_accessor(pos_handle.as<double>());

        const Eigen::Vector3d p0 = pos_acc.const_vector_attribute(vertices[0]);
        const Eigen::Vector3d p1 = pos_acc.const_vector_attribute(vertices[1]);
        const Eigen::Vector3d p2 = pos_acc.const_vector_attribute(vertices[2]);

        const Eigen::Vector3d p_mid = (p0 + p1 + p2) / 3.;

        const double area = (p1 - p0).cross(p2 - p0).norm();

        // Eigen::Vector3d n_face = ((p1 - p0).cross(p2 - p0)).normalized();
        //
        //// check face normal orientation
        //{
        //    Tuple tet_tuple = m.map_to_parent_tuple(s);
        //
        //    const auto tet_tag_acc =
        //    embedding_mesh.create_const_accessor<int64_t>(tet_tag_handle);
        //    // the tet tuple must be the one inside the offset (with a tag)
        //    if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
        //        tet_tuple = embedding_mesh.switch_tuple(tet_tuple, PrimitiveType::Tetrahedron);
        //        if (tet_tag_acc.const_scalar_attribute(tet_tuple) == 0) {
        //            log_and_throw_error("no tagged tet incident to offset");
        //        }
        //    }
        //    // get opposite vertex
        //    Tuple opposite_vertex = embedding_mesh.switch_tuples(
        //        tet_tuple,
        //        {PrimitiveType::Triangle, PrimitiveType::Edge, PrimitiveType::Vertex});
        //
        //    const auto tet_pos_acc =
        //        embedding_mesh.create_const_accessor<double>(embedding_pos_handle);
        //    const Eigen::Vector3d p3 = tet_pos_acc.const_vector_attribute(opposite_vertex);
        //
        //    const Eigen::Vector3d p_normal = p_mid + n_face;
        //    if (wmtk_orient3d(p0, p1, p2, p3) != wmtk_orient3d(p0, p1, p2, p_normal)) {
        //        n_face *= -1;
        //    }
        //    // n_face points now in the direction of the input
        //}

        constexpr double u = 0.1;
        const std::array<Eigen::Vector3d, 4> sample_points{
            {p_mid, (1 - u) * p0 + u * p_mid, (1 - u) * p1 + u * p_mid, (1 - u) * p2 + u * p_mid}};

        double sq_dist = -1;
        SimpleBVH::VectorMax3d nearest_point;

        auto get_point_and_normal =
            [&bvh, &offset_distance, &sq_dist, &nearest_point](
                const Eigen::Vector3d& _p) -> std::tuple<Eigen::Vector3d, Eigen::Vector3d> {
            bvh.nearest_facet(_p, nearest_point, sq_dist);

            const Eigen::Vector3d normal = (_p - nearest_point).normalized();
            const Eigen::Vector3d point = nearest_point + normal * offset_distance;
            return std::make_tuple(point, normal);
        };

        auto [q_point, q_normal] = get_point_and_normal(sample_points[0]);

        Quadrics q(q_point, q_normal);
        for (int64_t i = 1; i < sample_points.size(); ++i) {
            std::tie(q_point, q_normal) = get_point_and_normal(sample_points[i]);
            q += Quadrics(q_point, q_normal) * 0.1;
        }

        q *= area;

        m_matrix = q.m_matrix;

    } else if (s.primitive_type() == PrimitiveType::Edge) {
        assert(!m.is_multi_mesh_root());
        assert(!m.has_child_mesh());
        assert(m.top_simplex_type() == PrimitiveType::Edge);
        assert(pos_handle.dimension() == 2);

        const auto vertices = simplex::faces_single_dimension_tuples(m, s, PrimitiveType::Vertex);
        assert(vertices.size() == 2);

        const auto pos_acc = m.create_const_accessor(pos_handle.as<double>());

        const Eigen::Vector2d p0 = pos_acc.const_vector_attribute(vertices[0]);
        const Eigen::Vector2d p1 = pos_acc.const_vector_attribute(vertices[1]);

        const Eigen::Vector2d p_mid = (p0 + p1) / 2.;

        const double area = (p1 - p0).norm();

        constexpr double u = 0.1;
        const std::array<Eigen::Vector2d, 3> sample_points{
            {p_mid, (1 - u) * p0 + u * p_mid, (1 - u) * p1 + u * p_mid}};

        double sq_dist = -1;
        SimpleBVH::VectorMax3d nearest_point;

        auto get_point_and_normal =
            [&bvh, &offset_distance, &sq_dist, &nearest_point](
                const Eigen::Vector2d& _p) -> std::tuple<Eigen::Vector2d, Eigen::Vector2d> {
            const int64_t nearest_facet = bvh.nearest_facet(_p, nearest_point, sq_dist);

            const Eigen::Vector2d normal = (_p - nearest_point).normalized();
            const Eigen::Vector2d point = nearest_point + normal * offset_distance;
            return std::make_tuple(point, normal);
        };

        auto [q_point, q_normal] = get_point_and_normal(sample_points[0]);
        Eigen::Vector3d q_point_3(q_point[0], q_point[1], 0);
        Eigen::Vector3d q_normal_3(q_normal[0], q_normal[1], 0);

        Quadrics q(q_point_3, q_normal_3);
        for (int64_t i = 1; i < sample_points.size(); ++i) {
            std::tie(q_point, q_normal) = get_point_and_normal(sample_points[i]);
            Eigen::Vector3d q_point_3(q_point[0], q_point[1], 0);
            Eigen::Vector3d q_normal_3(q_normal[0], q_normal[1], 0);
            q += Quadrics(q_point_3, q_normal_3) * 0.1;
        }

        q *= area;

        m_matrix = q.m_matrix;
    } else if (s.primitive_type() == PrimitiveType::Vertex) {
        const auto cofaces = simplex::top_dimension_cofaces_tuples(m, s);
        for (const Tuple& f : cofaces) {
            const simplex::Simplex f_simplex(m, m.top_simplex_type(), f);
            (*this) += Quadrics(
                bvh,
                pos_handle,
                tet_tag_handle,
                bvh_mesh_position_handle,
                f_simplex,
                offset_distance);
        }
    } else {
        log_and_throw_error("Unknown primitive type for Quadrics.");
    }
}

Quadrics::Quadrics(const Eigen::Ref<const Eigen::Matrix4d> matrix)
    : m_matrix(matrix)
{}

Quadrics::Quadrics(const Quadrics& other)
    : m_matrix(other.m_matrix)
{}

Quadrics Quadrics::operator+=(const Quadrics& other)
{
    m_matrix += other.m_matrix;
    return *this;
}

Quadrics Quadrics::operator+(const Quadrics& other)
{
    Quadrics q(*this);
    q += other;
    return q;
}

Quadrics Quadrics::operator*(const double scalar)
{
    Quadrics q(*this);
    q.m_matrix *= scalar;
    return q;
}

Quadrics Quadrics::operator*=(const double scalar)
{
    m_matrix *= scalar;
    return *this;
}

Eigen::Vector3d Quadrics::solve(const Eigen::Ref<const Eigen::Vector3d> p)
{
    const Eigen::Matrix3d A = m_matrix.block(0, 0, 3, 3);
    const Eigen::Vector3d b = -m_matrix.block(0, 3, 3, 1);

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    svd.setThreshold(1e-2);

    // Lindstrom formula for QEM new position for singular matrices
    // from "Out-of-Corse Simplification of Large Polygonal Models"
    const Eigen::Vector3d v_svd = p + svd.solve(b - A * p);
    return v_svd;
}

double Quadrics::squared_distance_to_optimum(const Eigen::Ref<const Eigen::Vector3d> p)
{
    const Eigen::Vector3d p_opt = solve(p);
    return (p_opt - p).squaredNorm();
}

Eigen::VectorXd Quadrics::vectorized()
{
    return Eigen::Map<Eigen::VectorXd>(m_matrix.data(), 16);
}

double Quadrics::squared_distance_to_optimum(
    const attribute::MeshAttributeHandle& position_handle,
    const simplex::Simplex& s)
{
    const Mesh& m = position_handle.mesh();
    const auto acc = m.create_const_accessor(position_handle.as<double>());

    if (s.primitive_type() == PrimitiveType::Vertex) {
        return squared_distance_to_optimum(acc.const_vector_attribute(s.tuple()));
    }

    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    const auto vertices = simplex::faces_single_dimension_tuples(m, s, PrimitiveType::Vertex);
    for (const Tuple& v : vertices) {
        p += acc.const_vector_attribute(s.tuple());
    }
    p /= vertices.size();
    return squared_distance_to_optimum(p);
}

Eigen::Matrix4d Quadrics::attribute_to_matrix(
    const attribute::MeshAttributeHandle& quadric_handle,
    const Tuple& t)
{
    assert(quadric_handle.mesh().top_simplex_type() == PrimitiveType::Triangle);
    assert(quadric_handle.dimension() == 16);
    assert(quadric_handle.primitive_type() == PrimitiveType::Triangle);

    const Mesh& m = quadric_handle.mesh();
    const auto acc = m.create_const_accessor(quadric_handle.as<double>());
    const auto q = acc.const_vector_attribute(t);
    const Eigen::Matrix4d mat = Eigen::Map<const Eigen::Matrix4d>(q.data());

    return mat;
}

void Quadrics::matrix_to_attribute(
    attribute::MeshAttributeHandle& quadric_handle,
    const Tuple& t,
    const Eigen::Ref<Eigen::Matrix4d> mat)
{
    assert(quadric_handle.mesh().top_simplex_type() == PrimitiveType::Triangle);
    assert(quadric_handle.dimension() == 16);
    assert(quadric_handle.primitive_type() == PrimitiveType::Triangle);

    Mesh& m = quadric_handle.mesh();
    auto acc = m.create_accessor(quadric_handle.as<double>());
    acc.vector_attribute(t) = Eigen::Map<const Eigen::VectorXd>(mat.data(), 16);
}

double Quadrics::squared_distance_of_projected_simplex_to_optimum(
    const attribute::MeshAttributeHandle& position_handle,
    const attribute::MeshAttributeHandle& quadric_handle,
    const simplex::Simplex& s)
{
    Quadrics simplex_quadrics(quadric_handle, s);

    if (s.primitive_type() == PrimitiveType::Vertex) {
        return simplex_quadrics.squared_distance_to_optimum(position_handle, s);
    }

    const Mesh& m = position_handle.mesh();
    const auto acc = m.create_const_accessor(position_handle.as<double>());

    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    const auto vertices = simplex::faces_single_dimension(m, s, PrimitiveType::Vertex);
    for (const simplex::Simplex& v : vertices) {
        Quadrics vertex_quadrics(quadric_handle, v);
        p += vertex_quadrics.solve(acc.const_vector_attribute(v.tuple()));
    }
    p /= vertices.size();

    return simplex_quadrics.squared_distance_to_optimum(p);
}

} // namespace wmtk::components::internal::utils