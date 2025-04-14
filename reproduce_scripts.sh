mkdir -p output
cd output
echo ++++++++++ Puzzle ++++++++++
../build/applications/topological_offsets_app -j ../input_data/puzzle.json
echo ++++++++++ Hello World ++++++++++
../build/applications/topological_offsets_app -j ../input_data/hello_world.json
echo ++++++++++ Fertility Finite ++++++++++
../build/applications/topological_offsets_app -j ../input_data/fertility_finite.json
echo ++++++++++ Fertility Topo ++++++++++
../build/applications/topological_offsets_app -j ../input_data/fertility_topo.json