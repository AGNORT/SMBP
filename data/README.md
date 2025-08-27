## Description of Instances

Three sets of benchmark instances are used in this paper to test the proposed algorithms:

1. **Xu et al. (2023):** [[1]](#ref1) Instances for both SMKP and SMBP.
2. **Monaci et al. (2013):**[[2]](#ref2) Instances for SMKP.
3. **Ryu and Park (2021):** [[3]](#ref3)Additional instances for SMKP, generated following their method.

All instances share the same format, organized according to the following rules:

- **First line:**
   Contains three values:
  - The number of items (*n*),
  - The bin capacity,
  - The confidence level.
- **Next *n* lines:**
   Each line corresponds to one item and contains three parameters:
  1. Profit of the item (used only in SMKP),
  2. Average weight of the item,
  3. Weight variance of the item.

## Reference

<a id="ref1">[1]</a> Xu, L., d'Ambrosio, C., Haddad-Vanier, S., & Traversi, E. (2023). Branch and price for submodular bin packing. *EURO Journal on Computational Optimization*, *11*, 100074. (https://doi.org/10.1016/j.ejco.2023.100074).  

<a id="ref2">[2]</a> Monaci, M., Pferschy, U., & Serafini, P. (2013). Exact solution of the robust knapsack problem. *Computers & operations research*, *40*(11), 2625-2631. (https://doi.org/10.1016/j.cor.2013.05.005).  

<a id="ref3">[3]</a> Ryu, J., & Park, S. (2021). Robust solutions for stochastic and distributionally robust chance-constrained binary knapsack problems. *arXiv preprint arXiv:2105.11875*. ([ https://doi.org/10.48550/arXiv.2105.11875](https://doi.org/10.48550/arXiv.2105.11875)).  
