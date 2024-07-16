# `ChemicalBalancingUtilityNew`
A chemical balancing tool in C++ (new version) by **Frank Yang**

## A powerful chemical balancing tool
`ChemicalBalancingUtilityNew` (abbr. `CBU_Balancer`) is a C++ header that includes methods for balancing your chemical equations. With perfect encapsulation, you can efficiently and conveniently balance equations anywhere in your program. Our program also comes with a `CBU_Console` library, an original console designed for users who do not wish to delve deeply into C++.

## _Quick start_
The following program uses the `CBU_Console` class that comes with ,y program. The `boot()` method of this class will open a console where you can balance equations freely.
1. Create a new project and `main.cpp`. Be sure your project version is C++17 or later.
2. Download `CBU_Balancer.h` and `CBU_Console.h` and include them in your `main.cpp`:
   ```cpp
   #include "CBU_Balancer.h"
   #include "CBU_Console.h"
   ```
3. Call `CBU_Console` in your `main.cpp` and then you can use the powerful `ChemicalBalancingUtility`!
   ```cpp
   int main() {
     CBU_Console console = CBU_Console();
     console.boot();
     return 0;
   }
   ```
4. In the following console, type `quit()` to quit the console; type `multiple_results(off)` to disable multiple results, type `multiple_results(on)` to allow it. Multiple results is enabled by default (which is to provide at least one **absolutely correct answer**). **Directly type your chemical equation to call the built-in C`hemicalBalancingUtility` to balance**, for example,
   ```
   HNO3 -> NO2 + O2 + H2O
   ```


## Complete user guide (`CBU_Balancer`)
`CBU_Balancer` is the main class. This class allows you to balance equations that appear in your program (this is an extension of C++ functionality).

Use **`CBU_Balancer balancer = CBU_Balancer()`** to create a new `ChemicalBalancingUtilityNew` instance (all the needed features must be implemented through this instance);

### Public interfaces
1. `void balance(const std::string& equation)` receives a chemical equation as input and try to balance it. **The balancing data will be stored in `balancer`**. E.g.,
   ```cpp
   balancer.balance("C+O2->CO2");
   ```
2. `void balance_with_given_compounds(const std::vector<std::string>& reactants, const std::vector<std::string>& products)` two `std::vector<std::string>` objects (as reactants and products) and try to balance them. E.g.,
   ```cpp
   balancer.balance_with_given_compounds({"Zn", "HCl"}, {"ZnCl2", "H2"});
   ```
3. `std::pair<std::vector<std::string>, std::vector<std::string>> get_reactants_and_products()` returns an object where its `first` is the stored reactants vector and its `second` is the products vector. E.g.,
   ```cpp
   auto reactants_and_products = balancer.get_reactants_and_products();
   auto reactants = reactants_and_products.first;
   auto products = reactants_and_products.second;
   ```
4. `std::pair<std::vector<std::string>, std::vector<std::vector<unsigned>>> get_main_matrix()` returns the information of stored equation data. Its `first` is the elements occurred in the equation, and its `second` is the _main matrix_ of the equation. For the _main matrix_, its row is each element (sorted by `std::sort`) and its column is each compound (by given order);
5. `void clear_data()` clears all stored balancing data. **This is a must when you're to balance another equation.**.
6. The following shows an overall sample:
   ```
   balancer.balance("C+O2->CO2");
   std::cout << balancer.get_result() << std::endl;
   balancer.clear_data();
   balancer.balance_with_given_compounds({"Zn", "HCl"}, {"ZnCl2", "H2"});
   std::cout << balancer.get_main_matrix() << std::endl;
   std::cout << balancer.get_result() << std::endl;
   ```

## Configs (`CBU_Balancer`)
1. Use `set_multiple_results(bool option)` to set if you want a multiple results.
2. Use `set_max_coef(unsigned max_coef)` to set the maximum possible coefficient in the balanced equation. I suggest that `max_coef` should not be greater than 30.
3. Use `set_log_status(bool option)` to set if you want `std::clog` messages when have. 

