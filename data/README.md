## Description of Instances

Three sets of benchmark instances are used in this paper to test the proposed algorithms:

1. **Xu et al. (2023):** Instances for both SMKP and SMBP.
2. **Monaci et al. (2013):** Instances for SMKP.
3. **Ryu and Park (2021):** Additional instances for SMKP, generated following their method.

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
