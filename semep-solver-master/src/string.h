/**
 * Copyright (C) 2017 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gvpalma@usb.ve>
 */

#ifndef __STRING_H
#define __STRING_H

#define BUFSZ     256

struct char_array {
     unsigned nr;
     unsigned alloc;
     char *data;
};

struct string_array {
     unsigned int nr;
     unsigned int alloc;
     char **data;
};

static inline void add_char(struct char_array *buf, char ch)
{
     unsigned alloc, nr;

     alloc = buf->alloc;
     nr = buf->nr;
     if (nr == alloc) {
	  alloc = BUFSZ + alloc;
	  buf->data = realloc(buf->data, alloc);
	  if (!buf->data) {
	       fprintf(stderr, "Out of memory, realloc failed\n");
	       exit(1);
	  }
	  buf->alloc = alloc;
     }
     buf->data[nr] = ch;
     buf->nr++;
}

static inline void init_char_array(struct char_array *buf)
{
     buf->data = calloc(BUFSZ, 1);
     if (!buf->data) {
	  fprintf(stderr, "Out of memory, calloc failed\n");
	  exit(1);
     }
     buf->alloc = BUFSZ;
     buf->nr = 0;
}

static inline void string_clean(struct char_array *buf)
{
     buf->nr = 0;
}

static inline void free_string_array(struct string_array *sa)
{
     for (unsigned i = 0; i < sa->nr; i++)
	  if (sa->data[i])
	       free(sa->data[i]);
     free(sa->data);
     sa->nr = 0;
}

static inline void print_string_array(const struct string_array *sa)
{
     printf("\nNumber of elements %d\n", sa->nr);
     for (unsigned i = 0; i < sa->nr; i++)
	  printf("%s\n", sa->data[i]);
     printf("\n");
}

#endif /* __STRING_H */
