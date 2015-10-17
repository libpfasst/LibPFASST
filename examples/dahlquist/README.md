# Naming conventions

## *encap.f90*

- An encapsulation *encap* encapsulates data *dat*
- Data *dat* is comprised of
  - the number *nflds* of scalar fields;
  - the *shape* of the Cartesian grid, i.e.
    - the number *nx* of grid points in x-direction,
    - the number *ny* of grid points in y-direction, and
    - the number *nz* of grid points in z-direction; and
  - the *nflds*-by-*nx*-by-*ny*-by-*nz* array *arr* of scalar field values
- A *dat* qualified by, e.g., the color red, is denoted by *red_dat*
- A member such as *nflds* of a *dat* is referred to by *dat_nflds*
- A member such as *nflds* of a *dat* qualified by, e.g., the color red, is referred to by *red_dat_nflds*

## *transfer.f90*

- Input data defines running indices with respect to naming (e.g. for *nflds* the output data would we equally possible)
