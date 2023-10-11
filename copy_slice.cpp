#include <vector> 
#include <iostream>

int main(){
    std::vector<int> a{0,1,2,3,4,5,6,7,8,9};
    std::vector<int> b(10, 2);

    std::copy(a.begin() + 3, a.begin()+7, b.begin() + 6);
    for (int i: b) std::cout << i << ' ';
}