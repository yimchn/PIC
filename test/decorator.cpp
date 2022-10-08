#include <functional>
#include <iostream>

void test(int a, int b) { std::cout << a + b << std::endl; }

void do_previous_work() { std::cout << "do previous work\n"; }

void do_after_work() { std::cout << "do after work\n"; }

template <typename RET, typename... ARGS>
std::function<RET(ARGS... args)> decorator(RET (*p_func)(ARGS... args)) {
    return [=](ARGS... args) -> RET {
        do_previous_work();
        (*p_func)(args...);
        do_after_work();
    };
    //没有捕获任何变量的Lambda表达式可以转换成与它的调用原型一致的函数指针
}

int main() {
    int a = 10;
    int b = 20;
    std::function<void(int, int)> f;
    // f = test;
    // f(a, b);  // 第一次调用 f
    std::cout << "-----------------------------------------------" << std::endl;

    f = decorator(test);
    f(a, b);  // 第二次调用 f
    return 0;
}