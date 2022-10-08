#include <iostream>

template <typename T>
class Base {
   public:
    void interface() { static_cast<T*>(this)->imp(); }
    void imp() { std::cout << "in Base::imp" << std::endl; }
};

class Derived1 : public Base<Derived1> {
   public:
    void imp() { std::cout << "in Derived1::imp" << std::endl; }
};

class Derived2 : public Base<Derived2> {
   public:
    void imp() { std::cout << "in Derived2::imp" << std::endl; }
};

struct Derived3 : public Base<Derived3> {};

template <typename T>
void fun(T& base) {
    base.interface();
}

int main() {
    // Base<int> b;
    Derived1 d1;
    Derived2 d2;
    Derived3 d3;

    // fun(b);
    fun(d1);
    fun(d2);
    fun(d3);

    return 0;
}