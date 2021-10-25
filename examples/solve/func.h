
class Func {
    public:
        template<typename TypeVector>
        TypeVector operator()(const TypeVector& x, const double& t) const {
            return x; // Exponential derivative
        }
};

