/*
 * Copyright (c) 2016-2018 Ilya Kaliman
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>

#include "pt.h"

void *(*libpt_malloc)(size_t) = malloc;
void (*libpt_free)(void *) = free;

void
libpt_print_banner(void)
{
	printf("libpt (c) 2016-2018 Ilya Kaliman\n");
	printf("Fast Coupled Cluster Triples Corrections\n");
	printf("https://github.com/ilyak/libpt\n");
}

void
libpt_set_malloc(void *(*fn)(size_t))
{
	if (fn == NULL)
		fn = malloc;
	libpt_malloc = fn;
}

void
libpt_set_free(void (*fn)(void *))
{
	if (fn == NULL)
		fn = free;
	libpt_free = fn;
}
