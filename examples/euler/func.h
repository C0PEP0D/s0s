
template<typename TypeVector>
class Func {
    public:
        TypeVector operator()(const double* pX, const double& t) const {
            return TypeVector(pX); // Exponential derivative
        }
};

