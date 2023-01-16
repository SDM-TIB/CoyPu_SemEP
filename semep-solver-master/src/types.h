/**
 * Copyright (C) 2012-2017 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gvpalma@usb.ve>
 */

#ifndef ___TYPES_H
#define ___TYPES_H

struct int_array {
     unsigned int nr;
     unsigned int alloc;
     int *data;
};

struct matrix {
     unsigned int n;
     double **data;
};

#endif /* ___TYPES_H */
