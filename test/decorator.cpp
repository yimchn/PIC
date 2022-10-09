#include <functional>
#include <iostream>

/**
 * @brief 装饰器基类，内部封装一个函数器，用于指定要被包装的函数
 *
 * @tparam RET 被包装函数的返回值
 * @tparam ARGS 被包装函数的参数
 */
template <typename RET, typename... ARGS>
class Decorator {
   protected:
    std::function<RET(ARGS...)> decorated_func;  // 被包装的函数

   public:
    // 针对函数器的构造函数
    explicit Decorator(const std::function<RET(ARGS...)> &p_func)
        : decorated_func(p_func) {}

    // 针对函数指针的构造函数
    explicit Decorator(RET (*p_func)(ARGS...))
        : Decorator(std::function<RET(ARGS...)>(p_func)){};

    // 重载()运算符作为接口
    virtual RET operator()(ARGS... args) = 0;
};

#include <chrono>

/**
 * @brief 通过宏定义的匿名函数来处理函数返回值为void时编译报错的问题
 *
 * 当函数返回值为void时，模板替换的结果为void result =
 * this->decorated_func(args...); 由于C++中没有void类型变量，所以编译会报错。
 * 这里通过有返回值的lambda函数进行二次封装来解决这个问题，此处使用宏定义的形式。
 */
#define TIMING(content, text)   \
    Timer<void *, void *>(      \
        [&](void *) -> void * { \
            do {                \
                content         \
            } while (false);    \
            return nullptr;     \
        },                      \
        text)(nullptr)

template <typename RET, typename... ARGS>
class Timer : public Decorator<RET, ARGS...> {
   private:
    std::string tip;  // 针对计时器的说明字符串

   public:
    explicit Timer(const std::function<RET(ARGS...)> &p_func,
                   std::string tip = "")
        : Decorator<RET, ARGS...>(p_func), tip(std::move(tip)) {}

    explicit Timer(RET (*p_func)(ARGS...), std::string tip = "")
        : Decorator<RET, ARGS...>(p_func), tip(std::move(tip)) {}

    RET operator()(ARGS... args) override {
        // 记录函数运行前的时间
        auto start = std::chrono::steady_clock::now();

        // 调用被装饰的函数
        RET result = this->decorated_func(args...);

        // 记录函数运行后的时间
        auto end = std::chrono::steady_clock::now();
        // 计算时间差并输出
        auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << tip
                  << static_cast<double>(duration.count()) *
                         std::chrono::microseconds::period::num /
                         std::chrono::microseconds::period::den
                  << " s" << std::endl;

        return result;
    }
};

int add(int a, int b) {
    std::cout << a + b << std::endl;

    return a + b;
}

int main(int argc, char *argv[]) {
    // TIMING(add(1, 2);, "Test: ");
    Timer<void *, void *>(
        [&](void *) -> void * {
            do {
                add(1, 2);
            } while (false);
            return nullptr;
        },
        "Test: ")(nullptr);

    return 0;
};