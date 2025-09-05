#pragma once

#include <mshio/mshio.h>
#include <filesystem>
#include <wmtk/Mesh.hpp>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::utils {

class MshData
{
public:
    template <typename Fn>
    void add_tet_vertices(size_t num_vertices, const Fn& get_vertex_cb)
    {
        add_vertices<3>(num_vertices, get_vertex_cb);
    }

    template <typename Fn>
    void add_tets(size_t num_tets, const Fn& get_tet_cb)
    {
        add_simplex_elements<3>(num_tets, get_tet_cb);
    }

    template <int NUM_FIELDS, typename Fn>
    void add_tet_attribute(const std::string& name, const Fn& get_attribute_cb)
    {
        add_element_attribute<NUM_FIELDS, 3>(name, get_attribute_cb);
    }

    void save(const std::string& filename, bool binary = true)
    {
        m_spec.mesh_format.file_type = binary;
        mshio::validate_spec(m_spec);
        mshio::save_msh(filename, m_spec);
    }

private:
    template <int DIM, typename Fn>
    void add_vertices(size_t num_vertices, const Fn& get_vertex_cb)
    {
        static_assert(DIM >= 1 && DIM <= 3, "Only 1,2,3D elements are supported!");
        if (num_vertices == 0) return;
        mshio::NodeBlock block;
        block.num_nodes_in_block = num_vertices;
        block.tags.reserve(num_vertices);
        block.data.reserve(num_vertices * 3);
        block.entity_dim = DIM;
        block.entity_tag = m_spec.nodes.num_entity_blocks + 1;

        const size_t tag_offset = m_spec.nodes.max_node_tag;
        for (size_t i = 0; i < num_vertices; i++) {
            const auto& v = get_vertex_cb(i);
            block.tags.push_back(tag_offset + i + 1);
            block.data.push_back(v[0]);
            block.data.push_back(v[1]);
            block.data.push_back(v[2]);
        }

        m_spec.nodes.num_entity_blocks += 1;
        m_spec.nodes.num_nodes += num_vertices;
        m_spec.nodes.min_node_tag = 1;
        m_spec.nodes.max_node_tag += num_vertices;
        m_spec.nodes.entity_blocks.push_back(std::move(block));
    }

    template <int DIM, typename Fn>
    void add_simplex_elements(size_t num_elements, const Fn& get_element_cb)
    {
        static_assert(
            DIM == 1 || DIM == 2 || DIM == 3,
            "Only 1,2,3D simplex elements are supported");
        if (num_elements == 0) return;

        if (m_spec.nodes.num_nodes == 0) {
            throw std::runtime_error("Please add a vertex block before adding elements.");
        }
        const auto& vertex_block = m_spec.nodes.entity_blocks.back();
        assert(!vertex_block.tags.empty());
        if (vertex_block.entity_dim != DIM) {
            throw std::runtime_error("It seems the last added vertex block has different dimension "
                                     "than the elements you want to add.");
        }

        mshio::ElementBlock block;
        block.entity_dim = DIM;
        block.entity_tag = vertex_block.entity_tag;
        if constexpr (DIM == 1) {
            block.element_type = 1; // 2-node line.
        } else if constexpr (DIM == 2) {
            block.element_type = 2; // 3-node triangle.
        } else {
            block.element_type = 4; // 4-node tet.
        }
        block.num_elements_in_block = num_elements;

        const size_t vertex_offset = vertex_block.tags.front() - 1;
        const size_t tag_offset = m_spec.elements.max_element_tag;
        block.data.reserve(num_elements * (DIM + 2));
        for (size_t i = 0; i < num_elements; i++) {
            const auto& e = get_element_cb(i);
            block.data.push_back(tag_offset + i + 1); // element tag.
            for (size_t j = 0; j <= DIM; j++) {
                block.data.push_back(vertex_offset + e[j] + 1);
            }
        }

        m_spec.elements.num_entity_blocks++;
        m_spec.elements.num_elements += num_elements;
        m_spec.elements.min_element_tag = 1;
        m_spec.elements.max_element_tag += num_elements;
        m_spec.elements.entity_blocks.push_back(std::move(block));
    }

    template <int NUM_FIELDS, int ELEMENT_DIM, typename Fn>
    void add_element_attribute(const std::string& name, const Fn& get_attribute_cb)
    {
        static_assert(
            NUM_FIELDS == 1 || NUM_FIELDS == 3 || NUM_FIELDS == 9,
            "Only scalar, vector and tensor fields are supported as attribute!");
        if (m_spec.elements.entity_blocks.empty()) {
            throw std::runtime_error("Please add elements before adding element attributes!");
        }
        const auto& elem_block = m_spec.elements.entity_blocks.back();
        if (elem_block.entity_dim != ELEMENT_DIM) {
            throw std::runtime_error(
                "It seems the last added element block has different dimension "
                "than the element attribute you want to add.");
        }
        const size_t num_elements = elem_block.num_elements_in_block;

        mshio::Data data;
        data.header.string_tags = {name};
        data.header.real_tags = {0.0};
        data.header.int_tags = {0, NUM_FIELDS, int(num_elements), 0, ELEMENT_DIM};

        data.entries.resize(num_elements);
        for (size_t i = 0; i < num_elements; i++) {
            auto& entry = data.entries[i];
            entry.tag = elem_block.data[i * (ELEMENT_DIM + 2)];
            entry.data.reserve(NUM_FIELDS);
            const auto& attr = get_attribute_cb(i);
            if constexpr (NUM_FIELDS == 1) {
                entry.data.push_back(attr);
            } else {
                for (size_t j = 0; j < NUM_FIELDS; j++) {
                    entry.data.push_back(attr[j]);
                }
            }
        }

        m_spec.element_data.push_back(std::move(data));
    }

    template <int DIM>
    const mshio::NodeBlock* get_vertex_block() const
    {
        for (const auto& block : m_spec.nodes.entity_blocks) {
            if (block.entity_dim == DIM) {
                return &block;
            }
        }
        return nullptr;
    }

    template <int DIM>
    const mshio::ElementBlock* get_simplex_element_block() const
    {
        for (const auto& block : m_spec.elements.entity_blocks) {
            if (block.entity_dim == DIM) {
                return &block;
            }
        }
        return nullptr;
    }

private:
    mshio::MshSpec m_spec;
};

inline void write_msh(
    const attribute::MeshAttributeHandle& position_attribute,
    const attribute::MeshAttributeHandle& img_tag_attribute,
    const std::string& tag_attribute_name,
    const std::filesystem::path& file)
{
    const Mesh& m = position_attribute.mesh();

    if (m.top_simplex_type() != PrimitiveType::Tetrahedron) {
        log_and_throw_error("Can only write tet meshes to msh file.");
    }
    assert(position_attribute.primitive_type() == PrimitiveType::Vertex);
    assert(img_tag_attribute.primitive_type() == PrimitiveType::Tetrahedron);

    const auto pos_acc = m.create_const_accessor<double>(position_attribute);
    const auto tag_acc = m.create_const_accessor<int64_t>(img_tag_attribute);


    MshData msh;

    const auto& vtx = m.get_all(PrimitiveType::Vertex);
    msh.add_tet_vertices(vtx.size(), [&](size_t k) -> Eigen::VectorXd {
        return pos_acc.const_vector_attribute(vtx[k]);
    });

    const auto& tets = m.get_all(PrimitiveType::Tetrahedron);
    msh.add_tets(tets.size(), [&](size_t k) {
        const int64_t i = m.id(tets[k], PrimitiveType::Tetrahedron);
        const auto vs = m.orient_vertices(tets[k]);
        std::array<size_t, 4> data;
        for (int j = 0; j < 4; j++) {
            data[j] = m.id(vs[j], PrimitiveType::Vertex);
            assert(data[j] < vtx.size());
        }
        return data;
    });

    msh.add_tet_attribute<1>(tag_attribute_name, [&](size_t i) {
        return tag_acc.const_scalar_attribute(tets[i]);
    });

    msh.save(file.string(), true);
}

} // namespace wmtk::utils