//
// Created by nack2 on 09.06.2025.
//

#include "../library.h" //do not copy, link in a standard way
#include <chrono>
#include <iostream>

///Code to copy to your main
struct fun
{
    std::function<double(double )> f;
    std::pair<double,double> intervals;
};

fun f1([](double x){return pow(x,2)*pow(sin(x),3); },std::make_pair(1,4.764798248));
fun f2([](double x){return exp(pow(x,2))*(1-x); },std::make_pair(-2,3.20870913294));

int main() {
    double cal1naswer = -10.1010101105917;
    double cal2naswer = -9876.54321007546;
    double cal3answer = 4.2025;


    for (int i =10; i < 101 ; i=i+10)
    {
        int precison = i;
        std::cout <<i << std::endl;
        numericIntegration n1(4,f1.intervals,f1.f);
        auto t1 = std::chrono::high_resolution_clock::now();
        double cal1 = n1.gaussLegendre(precison);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto ms_int = duration_cast<std::chrono::nanoseconds>(t2 - t1);
        std::cout<<"Wynik dla calki 1: "<<cal1<<" w czasie: "<<ms_int<<" milisekund\n";
        std::cout<<"Blad:"<<sqrt(pow(cal1 - cal1naswer,2))<<std::endl;

        t1 = std::chrono::high_resolution_clock::now();
        numericIntegration n2(4,f2.intervals,f2.f);
        double cal2 = n2.gaussLegendre(precison);
        t2 = std::chrono::high_resolution_clock::now();
        ms_int = duration_cast<std::chrono::nanoseconds>(t2 - t1);
        std::cout<<"Wynik dla calk 2: "<<cal2<<" w czasie: "<<ms_int<<" milisekund\n";
        std::cout<<"Blad:"<<sqrt(pow(cal2 - cal2naswer,2))<<std::endl;
    }

    return 0;
}
