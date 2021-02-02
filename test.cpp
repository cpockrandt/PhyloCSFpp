#include <random>
#include <iostream>

int main()
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(1.0, 2.0);
    for (size_t i = 0; i < 10; ++i)
        std::cout << dis(gen) << '\t';
    std::cout << '\n';
}
