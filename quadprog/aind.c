/*  Copyright (C) 1997-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com> */

/*  This program is free software; you can redistribute it and/or modify */
/*  it under the terms of the GNU General Public License as published by */
/*  the Free Software Foundation; either version 2 of the License, or */
/*  (at your option) any later version. */

/*  This program is distributed in the hope that it will be useful, */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*  GNU General Public License for more details. */

/*  You should have received a copy of the GNU General Public License */
/*  along with this program; if not, write to the Free Software */
/*  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, */
/*  USA. */

/*  this routine checks whether Aind has valid entries, i.e., */
/*    1) 1<= Aind(1,i) <= n for i=1,...,q (number of constraints) */
/*    2) 1<= Aind(j,i) <= n for j=2,...,Aind(1,i)+1, i=1,...,q */

/*  Aind is a m times q matrix constructed in Splus */

/* Subroutine */ int aind_(int *ind, int *m, int *q, int *n,
	int *ok)
{
    /* System generated locals */
    int ind_dim1, ind_offset, i__1, i__2;

    /* Local variables */
    int i__, j;

    /* Parameter adjustments */
    ind_dim1 = *m;
    ind_offset = 1 + ind_dim1;
    ind -= ind_offset;

    /* Function Body */
    *ok = 0;
    i__1 = *q;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ind[i__ * ind_dim1 + 1] < 1 || ind[i__ * ind_dim1 + 1] > *n) {
	    return 0;
	}
	i__2 = ind[i__ * ind_dim1 + 1] + 1;
	for (j = 2; j <= i__2; ++j) {
	    if (ind[j + i__ * ind_dim1] < 1 || ind[j + i__ * ind_dim1] > *n) {
		return 0;
	    }
	}
    }
    *ok = 1;
    return 0;
} /* aind_ */

