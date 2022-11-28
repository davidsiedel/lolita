auto constexpr cell_field = L2("CellDisplacement", 2, Monomial(1));
auto constexpr face_field = L2("FaceDisplacement", 1, Monomial(1));
auto constexpr cg_field = H1("Displacemnent", Lagrange(1));
auto constexpr dg_field = L2("Displacemnent", Lagrange(2));

Gradient(cell_field);
Identity(cell_field);
Trace(cell_field, 2);