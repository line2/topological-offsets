#pragma once
#include <nlohmann/json.hpp>
namespace {

nlohmann::json topological_offsets_spec = R"(
[
  {
    "pointer": "/",
    "type": "object",
    "required": ["input", "distance"],
    "optional": [
      "2d",
      "output",
      "tag_input",
      "tag_attribute",
      "finite_offset",
      "passes",
      "embedding_passes",
      "offset_passes",
      "min_edge_length",
      "max_edge_length",
      "min_normal_deviation",
      "max_normal_deviation",
      "side",
      "convergence_max",
      "distance_function",
      "relative_distance_and_length",
      "intermediate_output",
      "use_simplicial_embedding",
      "use_warm_start",
      "expand_topo",
      "adapt_offset_distance",
      "restrict_min_edge_length_to_input_avg",
      "invert_tags",
      "report",
      "log_file",
      "input_path",
      "opt_log_level",
      "DEBUG_print_embedding",
      "DEBUG_print_offset",
      "DEBUG_print_smooth"
    ]
  },
  {
    "pointer": "/input",
    "type": "string"
  },
  {
    "pointer": "/distance",
    "type": "list",
    "min": 1,
    "doc": "List of desired offset distances."
  },
  {
    "pointer": "/distance/*",
    "type": "float",
    "min": 0,
    "doc": "Desired offset distance."
  },
  {
    "pointer": "/2d",
    "type": "bool",
    "default": false,
    "doc": "The input is 2d."
  },
  {
    "pointer": "/output",
    "type": "string",
    "default": ""
  },
  {
    "pointer": "/tag_input",
    "type": "string",
    "default": ""
  },
  {
    "pointer": "/tag_attribute",
    "type": "string",
    "default": ""
  },
  {
    "pointer": "/finite_offset",
    "type": "bool",
    "default": false,
    "doc": "If true, the component will compute a finite offset instead of the topological one."
  },
  {
    "pointer": "/min_edge_length",
    "type": "float",
    "default": null,
    "doc": "Edges should not be shorter than this value. If chosen too large it might conflict with the normal deviation. This value is optional. If empty, it is computed based on offset distance and max normal deviation."
  },
  {
    "pointer": "/max_edge_length",
    "type": "float",
    "default": null,
    "doc": "Edges should not be larger than this value. This value is optional. If empty, it is set to infinity."
  },
  {
    "pointer": "/min_normal_deviation",
    "type": "float",
    "default": 2,
    "doc": "The optimization will increase the local target edge length if triangles have a normal deviation below the minimum."
  },
  {
    "pointer": "/max_normal_deviation",
    "type": "float",
    "default": 15,
    "doc": "The optimization will decrease the local target edge length if triangles have a normal deviation above the maximum."
  },
  {
    "pointer": "/side",
    "type": "string",
    "options": ["in", "out", "double"],
    "default": "out",
    "doc": "On which side of the input the offset should be generated: in = inside, out = outside, double = both sides."
  },
  {
    "pointer": "/convergence_max",
    "type": "float",
    "default": 0.5,
    "doc": "The maximum value difference before convergence."
  },
  {
    "pointer": "/distance_function",
    "type": "string",
    "default": "",
    "doc": "It is possible to initialize the sizing field by defining a volumetric distance function, e.g., x*y*z. Available variables: x,y,z, and d as distance."
  },
  {
    "pointer": "/passes",
    "type": "int",
    "default": 10
  },
  {
    "pointer": "/embedding_passes",
    "type": "int",
    "default": 2
  },
  {
    "pointer": "/offset_passes",
    "type": "int",
    "default": 2
  },
  {
    "pointer": "/relative_distance_and_length",
    "type": "bool",
    "default": true
  },
  {
    "pointer": "/intermediate_output",
    "type": "bool",
    "default": false
  },
  {
    "pointer": "/use_simplicial_embedding",
    "type": "bool",
    "default": true
  },
  {
    "pointer": "/use_warm_start",
    "type": "bool",
    "default": false,
    "doc": "Add all tets to the topological offset input that are within offset distance. This option has no effect when a finite offset is generated."
  },
  {
    "pointer": "/expand_topo",
    "type": "bool",
    "default": true,
    "doc": "Expand the topological offset by adding tets without changing the offsets' topology. This option has no effect when a finite offset is generated."
  },
  {
    "pointer": "/adapt_offset_distance",
    "type": "bool",
    "default": true,
    "doc": "In case of the warm start, locally adapt the offset distance to avoid colliding and potentially degenerating offsets."
  },
  {
    "pointer": "/restrict_min_edge_length_to_input_avg",
    "type": "bool",
    "default": false,
    "doc": "The minimal edge length will be bounded by the average edge length of the input. This setting increases the error in the normal deviation but reduces runtime."
  },
  {
    "pointer": "/invert_tags",
    "type": "bool",
    "default": false,
    "doc": "Switch inside and outside tags on top dimension simplices."
  },
  {
    "pointer": "/report",
    "type": "string",
    "default": "topological_offsets_report.json"
  },
  {
    "pointer": "/log_file",
    "type": "string",
    "default": "topological_offsets_log.txt"
  },
  {
    "pointer": "/input_path",
    "type": "string",
    "default": ""
  },
  {
    "pointer": "/opt_log_level",
    "type": "int",
    "default": 6
  },
  {
    "pointer": "/DEBUG_print_embedding",
    "type": "bool",
    "default": false,
    "doc": "Print intermediate embedding results"
  },
  {
    "pointer": "/DEBUG_print_offset",
    "type": "bool",
    "default": false,
    "doc": "Print intermediate offset results"
  },
  {
    "pointer": "/DEBUG_print_smooth",
    "type": "bool",
    "default": false,
    "doc": "Print intermediate smooth results"
  }
]
)"_json;

}
