# `ChemicalBalancingUtilityNew`
A chemical balancing tool in C++ (new version) by **Frank Yang**

## A powerful chemical balancing tool
`ChemicalBalancingUtilityNew` is a C++ library that includes methods for balancing your chemical equations. With perfect encapsulation, you can efficiently and conveniently balance equations anywhere in your program. Our program also comes with a `Console` library, an original console designed for chemists who do not wish to delve deeply into C++. 

## How to use `ChemicalBalancingUtilityNew` library
1. Use **`CBU_Balancer balancer = CBU_Balancer()`** to create a new `ChemicalBalancingUtilityNew` instance (all the needed features must be implemented through this instance);
2. Now you prepare a `std::string` object which is your chemical equation **(the chemical equation must be like `CaO + H2O -> Ca(OH)2`, any other form except the existance of the spacebar is not allowed)**, we assume an object `std::string eq = "CaO + H2O -> Ca(OH)2"`;
3. **Use `balancer.balance(eq)` to balance your equation.** Now the result is stored in the `balancer`.
4. **Use `balancer.get_result()` to get out the result (in `std::string` form).** As the example of `eq`, the printed result should be `"CaO + H2O == Ca(OH)2"`.

