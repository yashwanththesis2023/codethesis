#include <iostream>
#include <stdexcept>

int divide(int dividend, int divisor) {
     if (divisor == 0) {
         throw std::string("Division by zero");
     }

     return dividend / divisor;
}

int main() {
     try {
         int result = divide(10, 0);
         std::cout << "Result: " << result << std::endl;
     } catch (const std::string& error) {
         std::cout << "An exception occurred: " << error << std::endl;
     }

     return 0;
}