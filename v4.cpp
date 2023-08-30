#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    int row, column;
    double data[];
} matrix;

struct matrix_list {
	matrix *mat;
	matrix_list *next;
};

typedef matrix_list *ctrlpoint_list;

matrix *matrix_create(int row, int column);
matrix *matrix_create_identity(int n);
matrix *matrix_clone(matrix *mat);
void matrix_del(matrix *mat);
double *matrix_at(matrix *mat,int i, int j);
matrix *matrix_invert(matrix *mat);
matrix *matrix_multi(matrix *a, matrix *b);
matrix *matrix_pinv(matrix *mat);
void matrix_fprint(matrix *mat, FILE *stream);
void points_fprint(char name, int start, matrix* mat, FILE *stream);
void matrix_fscan(matrix *mat, FILE *stream);
matrix *matrix_bernstein(matrix *t, int n);
matrix *chordlength(matrix *d, int first, int last);
double distance(matrix *d, int a, int b);
double combination(double k, double n);
matrix *matrix_sub(matrix *a, matrix *b);
matrix *matrix_delta(int n);
double max(double a, double b);
double norm(matrix *mat);
void gauss_newton(matrix *r, matrix *r_deriv, matrix *t);
matrix *fitbezier(matrix *d, int n);
matrix *generatecubic(matrix *d, int first, int last, matrix *b, matrix *v);
void fitcubic(matrix *d, int first, int last, matrix *v, ctrlpoint_list &ctrlp_list);
matrix *matrix_r_maxerror(matrix *b, matrix *p, matrix *d, int first, int last, double &maxerror, int &splitpoint);
ctrlpoint_list fitcurve_cubic(matrix *d);
void addlist(matrix_list *&list, matrix* x);
void dellist(matrix_list *&list);
void input(matrix *&d, int &which_method);
void getdata(matrix *&d, FILE *stream);
FILE *letsopen(int which_file);
void output(matrix_list *list, int which_method, FILE *stream);
double error = 0.2, eps = 1e-3, interationerror = 5, maxinterations = 100;

int main() {
	printf("*****XAP XI DUONG CONG BEZIER BANG BINH PHUONG TOI THIEU*****\n");
	printf("_____________________________________________________________\n\n");
	char quit = 'c';
	int n, which_method;
	matrix *d, *p;
	ctrlpoint_list ctrlp_list = NULL;
	FILE *f;
	while(quit != 'q') {
		input(d, which_method);
		f = fopen("output.txt", "w");
		if(which_method == 1) {
			printf("Nhap bac cua duong cong(n < %d): ", d -> row);
			scanf("%d", &n);
			if(n < d -> row){
				p = fitbezier(d, n);
				addlist(ctrlp_list, p);
				output(ctrlp_list, which_method, stdout);
				output(ctrlp_list, which_method, f);
				dellist(ctrlp_list);				
			}
		}
		else {
			ctrlp_list = fitcurve_cubic(d);
			output(ctrlp_list, which_method, stdout);
			output(ctrlp_list, which_method, f);
			dellist(ctrlp_list);
		}
		fclose(f);
		matrix_del(d);
		printf("Nhap [q] de thoat: ");
		fflush(stdin);
		scanf("%c", &quit);
	}
	return 0;
}

// Tao mot ma tran rong
matrix *matrix_create(int row, int column) {
    matrix *mat = (matrix*)malloc(sizeof(matrix) + sizeof(double) * row * column);
    mat -> row = row;
    mat -> column = column;
    return mat;
}

// Tao ma tran don vi
matrix *matrix_create_identity(int n) {
    matrix *mat = matrix_create(n, n);
    memset(mat -> data, 0, sizeof(double) * n * n);
    for(int i = 0; i < n; i++) *matrix_at(mat, i, i) = 1;
    return mat;
}

// Copy ma tran
matrix *matrix_clone(matrix *mat) {
    matrix *clone = matrix_create(mat -> row, mat -> column);
    memcpy(clone -> data, mat -> data, sizeof(double) * mat -> row * mat -> column);
    return clone;
}

// Xoa ma tran(giai phong bo nho)
void matrix_del(matrix *mat) {
    free(mat);
}

// Tra ve dia chi mot phan tu cua ma tran
double *matrix_at(matrix *mat, int i, int j) {
    return &(mat -> data[i * mat -> column + j]);
}

// Nghich dao ma tran su dung phuong phap khu Gauss
matrix *matrix_invert(matrix *input) {
    if(input -> row != input -> column)
    	return NULL;
	matrix *mat = matrix_clone(input), *inv = matrix_create_identity(input -> row), *tmp = matrix_create(1, input -> column);
    for(int i = 0; i < input -> column; i++) {
        if(*matrix_at(mat, i, i) == 0) {
            int swap = i;
            for(; *matrix_at(mat, i, swap) == 0; swap++);
            if(swap == input -> column) {
                matrix_del(mat);
                matrix_del(tmp);
                matrix_del(inv);
                return NULL;
            }
            memcpy(matrix_at(tmp, 0, 0),    matrix_at(mat, 0, i),    sizeof(double) * input -> column);
            memcpy(matrix_at(mat, 0, i),    matrix_at(mat, 0, swap), sizeof(double) * input -> column);
            memcpy(matrix_at(mat, 0, swap), matrix_at(tmp, 0, 0),    sizeof(double) * input -> column);
            memcpy(matrix_at(tmp, 0, 0),    matrix_at(inv, 0, i),    sizeof(double) * input -> column);
            memcpy(matrix_at(inv, 0, i),    matrix_at(inv, 0, swap), sizeof(double) * input -> column);
            memcpy(matrix_at(inv, 0, swap), matrix_at(tmp, 0, 0),    sizeof(double) * input -> column);
        }
        for(int x = 0; x < input -> column; x++) {
            if(x > i)
				*matrix_at(mat, x, i) /= *matrix_at(mat, i, i);
            *matrix_at(inv, x, i) /= *matrix_at(mat, i, i);
        }
        *matrix_at(mat, i, i) = 1;
        for(int y = 0; y < input -> column; y++) {
            if(y == i) continue;
            if(*matrix_at(mat, i, y) == 0) continue;
            double mul = - *matrix_at(mat, i, y);
            for(int x = 0; x < input -> column; x++) {
                if(x >= i) *matrix_at(mat, x, y) += mul * *matrix_at(mat, x, i);
                *matrix_at(inv, x, y) += mul * *matrix_at(inv, x, i);
            }
        }
    }
    matrix_del(tmp);
    matrix_del(mat);
    return inv;
}

// Nhan ma tran
matrix *matrix_multi(matrix *a, matrix *b) {
	if(a -> column != b -> row)
		return NULL;
	matrix *tmp = matrix_create(a -> row, b -> column);
	memset(tmp -> data, 0, sizeof(double) * tmp -> row * tmp -> column);
	for(int i1 = 0; i1 < a -> row; i1++)
		for(int j2 = 0; j2 < b -> column; j2++)
			for(int j1 = 0; j1 < a -> column; j1++)
				*matrix_at(tmp, i1, j2) += (*matrix_at(a, i1, j1))*(*matrix_at(b, j1, j2));
	return tmp;
}

// Nghich dao gia(pseudo inverse)
matrix *matrix_pinv(matrix *mat) {
	matrix *mtm = matrix_create(mat -> column, mat -> column);
	memset(mtm -> data, 0, sizeof(double) * mtm -> row * mtm -> column);
	for(int i = 0; i < mtm -> row; i++)
		for(int j = 0; j < mtm -> column; j++)
			for(int k = 0; k < mat -> row; k++)
				*matrix_at(mtm, i, j) += *matrix_at(mat, k, i) * *matrix_at(mat, k, j);
	matrix *inv = matrix_invert(mtm); matrix_del(mtm);
	if(inv == NULL)
		return NULL;
	matrix *pinv = matrix_create(inv -> row, mat -> row);
	memset(pinv -> data, 0, sizeof(double) * pinv -> row * pinv -> column);
	for(int i = 0; i < pinv -> row; i++)
		for(int j = 0; j < pinv -> column; j++)
			for(int k = 0; k < mat -> column; k++)
				*matrix_at(pinv, i, j) += *matrix_at(inv, i, k) * *matrix_at(mat, j, k);
	return pinv;
}

// In ma tran
void matrix_fprint(matrix *mat, FILE *stream) {
	for(int i = 0; i < mat -> row; i++) {
		for(int j = 0; j < mat -> column; j++)
			fprintf(stream, "%+6lf\t", *matrix_at(mat, i, j));
		fprintf(stream, "\n");
    }
}

// In cac diem
void points_fprint(char name, int start, matrix* mat, FILE *stream) {
	for(int i = 0; i < mat -> row; i++) {
		fprintf(stream, "%c[%d]: [", name, start + i);
		for(int j = 0; j < mat -> column; j++) {
			fprintf(stream, "%+6lf", *matrix_at(mat, i, j));
			if(j != mat -> column - 1)
				fprintf(stream, ", ");
		}
		fprintf(stream, "]\n");
    }
}

// Nhap ma tran
void matrix_fscan(matrix *mat, FILE *stream) {
	for(int i = 0; i < mat -> row; i++)
		for(int j = 0; j < mat -> column; j++)
			fscanf(stream, "%lf", matrix_at(mat, i, j));
}

// To hop chap k cua n
double combination(double k, double n) {
	if(k > n) return -1;
	if(k == 0 || k == n) return 1;
	if(k == 1 || k == n-1) return n;
	return combination(k, n-1) * n / (n-k);
}

// Tao ma tran bernstein
matrix *matrix_bernstein(matrix *t, int n) {
	matrix *t_ = matrix_create(t -> row, n + 1);
	for(int i = 0; i < t -> row; i++) {
		*matrix_at(t_, i, n) = 1;
		for(int j = n-1; j >= 0; j--)
			*matrix_at(t_, i, j) = (*matrix_at(t_, i, j+1)) * (*matrix_at(t, i, 0));
	}
	matrix *m = matrix_create(n + 1, n + 1);
	memset(m -> data, 0, sizeof(double) * (n + 1) * (n + 1));
	for(int i = 0; i <= n; i++) {
		double x = combination(i, n);
		for(int j = 0; j <= n - i; j++)
			*matrix_at(m, i, j) = x * combination(n - i - j, n - i) * (((n - i - j) % 2 == 0)?1:-1);
	}
	matrix *b = matrix_multi(t_, m);
	matrix_del(t_);
	matrix_del(m);
	return b;
}

// Khoang cach giua 2 diem
double distance(matrix *d, int a, int b) {
	double dist = 0;
	for(int j = 0; j < d -> column; j++)
		dist += (*matrix_at(d, b, j) - *matrix_at(d, a, j)) * (*matrix_at(d, b, j) - *matrix_at(d, a, j));
	return sqrt(dist);
}

// Tao t bang phuong phap day cung(chord length)
matrix *chordlength(matrix *d, int first, int last) {
	double total = 0;
	int m = last - first + 1;
	for(int i = 1; i < m; i++)
		total += distance(d, i - 1, i);
	matrix *u = matrix_create(m, 1);
	*matrix_at(u, 0, 0) = 0;
	for(int i = 1; i < m; i++)
		*matrix_at(u, i, 0) = *matrix_at(u, i - 1, 0) + distance(d, i - 1, i) / total; 
	return u;
}

// Tru 2 ma tran
matrix *matrix_sub(matrix *a, matrix *b) {
	if(a -> row != b -> row || a -> column != b -> column)
		return NULL;
	matrix *tmp = matrix_create(a -> row, a -> column);
	for(int i = 0; i < a -> row; i++)
		for(int j = 0; j < a -> column; j++)
			*matrix_at(tmp, i, j) = *matrix_at(a, i, j) - *matrix_at(b, i, j);
	return tmp;
}

//Tao ma tran delta tinh dao ham
matrix *matrix_delta(int n) {
	matrix *delta = matrix_create(n, n + 1);
	memset(delta -> data, 0, sizeof(double) * n * (n + 1));
	for(int i = 0; i < n; i++) {
		*matrix_at(delta, i, i) = -n;
		*matrix_at(delta, i, i + 1) = n;
	}
	return delta;
}

// Tinh R'
matrix *r_derivative(matrix *p, matrix *t, int n) {
	matrix *delta, *tmp_1, *tmp_2, *r_deriv;
	delta = matrix_delta(n);
	tmp_1 = matrix_bernstein(t, n - 1);
	tmp_2 = matrix_multi(tmp_1, delta);
	r_deriv = matrix_multi(tmp_2, p);
	matrix_del(delta);
	matrix_del(tmp_1);
	matrix_del(tmp_2);
	return r_deriv;
}

// Cai thien t bang phuong phap Gauss - Newton
void gauss_newton(matrix *r, matrix *r_deriv, matrix *t) {
	double tmp_1, tmp_2, min, max;
	for(int i = 0; i < r -> row; i++) {
		tmp_1 = tmp_2 = 0;
		for(int j = 0; j < r -> column; j++) {
			tmp_1 += *matrix_at(r_deriv, i, j) * *matrix_at(r_deriv, i, j);
			tmp_2 += *matrix_at(r_deriv, i, j) * *matrix_at(r, i, j);
		}
		*matrix_at(t, i, 0) -= tmp_2 / tmp_1;
	}
	min = max = *matrix_at(t, 0 ,0);
	for(int i = 1; i < t -> row; i++) {
		if(*matrix_at(t, i, 0) < min)
			min = *matrix_at(t, i, 0);
		if(*matrix_at(t, i, 0) > max)
			max = *matrix_at(t, i, 0);
	}
	max -= min;
	for(int i = 0; i < t -> row; i++) {
		*matrix_at(t, i, 0) -= min;
		*matrix_at(t, i, 0) /= max;
	}
}

// Frobenius norm
double norm(matrix *mat) {
	double tmp = 0;
	for(int i = 0; i < mat -> row; i++)
		for(int j = 0; j < mat -> column; j++)
			tmp += (*matrix_at(mat, i, j)) * (*matrix_at(mat, i, j));
	return sqrt(tmp);
}

// Tra ve gia tri max
double max(double a, double b) {
	if(a > b) return a;
	return b;
}

// Ham fit voi bac n chi dinh
matrix *fitbezier(matrix *d, int n) {
	matrix *t, *b, *pinv_b, *p, *tmp, *r_new, *r_old, *r_deriv, *stop; 
	t = chordlength(d, 0, d -> row - 1);
	b = matrix_bernstein(t, n);
	pinv_b = matrix_pinv(b);
	p = matrix_multi(pinv_b, d); matrix_del(pinv_b);
	
	tmp = matrix_multi(b, p); matrix_del(b);
	r_new = matrix_sub(tmp, d); matrix_del(tmp);
	r_old = matrix_create(r_new -> row, r_new -> column);
	memset(r_old -> data, 0, sizeof(double) * r_old -> row * r_old -> column);
	stop = matrix_sub(r_new, r_old);
	while(norm(stop) / max(norm(r_new), 1) > eps) {
		r_deriv = r_derivative(p, t, n); matrix_del(p);
		gauss_newton(r_new, r_deriv, t); matrix_del(r_deriv);

		b = matrix_bernstein(t, n);
		pinv_b = matrix_pinv(b);
		p = matrix_multi(pinv_b, d); matrix_del(pinv_b);
		
		matrix_del(r_old);
		r_old = matrix_clone(r_new); matrix_del(r_new);
		tmp = matrix_multi(b, p); matrix_del(b);
		r_new = matrix_sub(tmp, d); matrix_del(tmp);
		matrix_del(stop);
		stop = matrix_sub(r_new, r_old);
	}
	matrix_del(stop); matrix_del(r_old); matrix_del(t); matrix_del(r_new);
	return p;
}

// Ham fit bac 3 voi 2 vector
matrix *generatecubic(matrix *d, int first, int last, matrix *b, matrix *v) {
	int m = last - first + 1;
	matrix *p = matrix_create(4, d -> column);
	for(int j = 0; j < p -> column; j++) {
		*matrix_at(p, 0, j) = *matrix_at(d, first, j);
		*matrix_at(p, 1, j) = *matrix_at(d, first, j);
		*matrix_at(p, 2, j) = *matrix_at(d, last, j);
		*matrix_at(p, 3, j) = *matrix_at(d, last, j);
	}
	matrix *a_1 = matrix_create(m, d -> column);
	matrix *a_2 = matrix_create(m, d -> column);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < d -> column; j++) {
			*matrix_at(a_1, i, j) = *matrix_at(v, 0, j) * *matrix_at(b, i, 1);
			*matrix_at(a_2, i, j) = *matrix_at(v, 1, j) * *matrix_at(b, i, 2);
		}
	matrix *c = matrix_create(2, 2); memset(c -> data, 0, sizeof(double) * 4);
	matrix *x = matrix_create(2, 1); memset(x -> data, 0, sizeof(double) * 2);
	double tmp;
	for(int i = 0; i < m; i++)
		for(int j = 0; j < d -> column; j++) {
			*matrix_at(c, 0, 0) += *matrix_at(a_1, i, j) * *matrix_at(a_1, i, j);
			*matrix_at(c, 1, 1) += *matrix_at(a_2, i, j) * *matrix_at(a_2, i, j);
			*matrix_at(c, 0, 1) = *matrix_at(c, 1, 0) += *matrix_at(a_1, i, j) * *matrix_at(a_2, i, j);
			tmp = *matrix_at(d, first + i, j) - *matrix_at(d, first, j) * (*matrix_at(b, i, 0) + *matrix_at(b, i, 1)) - *matrix_at(d, last, j) * (*matrix_at(b, i, 2) + *matrix_at(b, i, 3));
			*matrix_at(x, 0, 0) += tmp * *matrix_at(a_1, i, j);
			*matrix_at(x, 1, 0) += tmp * *matrix_at(a_2, i, j);
		}
	matrix *inv_c = matrix_invert(c);
	matrix *alpha;
	if(inv_c == NULL) {
		alpha = matrix_create(2, 1);
		memset(alpha -> data, 0, 2 * sizeof(double));
	}
	else
		alpha = matrix_multi(inv_c, x);
	double seglength = distance(d, first, last);
	double epsilon = seglength * 1.0e-6;
	if(*matrix_at(alpha, 0, 0) < epsilon || *matrix_at(alpha, 1, 0) < epsilon) {
		double dist = seglength / 3.0;
		for(int j = 0; j < p -> column; j++) {
			*matrix_at(p, 1, j) += *matrix_at(v, 0, j) * dist;
			*matrix_at(p, 2, j) += *matrix_at(v, 1, j) * dist;
		}
		matrix_del(a_1); matrix_del(a_2); matrix_del(c); matrix_del(x); matrix_del(inv_c); matrix_del(alpha);
		return p;
	}
	for(int j = 0; j < p -> column; j++) {
		*matrix_at(p, 1, j) += *matrix_at(alpha, 0, 0) * *matrix_at(v, 0, j);
		*matrix_at(p, 2, j) += *matrix_at(alpha, 1, 0) * *matrix_at(v, 1, j);
	}
	matrix_del(a_1); matrix_del(a_2); matrix_del(c); matrix_del(x); matrix_del(inv_c); matrix_del(alpha);
	return p;
}

// Ham fit tap diem du lieu voi nhieu duong cong bezier bac 3
void fitcubic(matrix *d, int first, int last, matrix *v, ctrlpoint_list &ctrlp_list) {
	if(last - first == 1) {
		matrix *p = matrix_create(4, d -> column);
		for(int j = 0; j < p -> column; j++) {
			*matrix_at(p, 0, j) = *matrix_at(d, first, j);
			*matrix_at(p, 1, j) = *matrix_at(d, first, j);
			*matrix_at(p, 2, j) = *matrix_at(d, last, j);
			*matrix_at(p, 3, j) = *matrix_at(d, last, j);
		}
		double dist = distance(d, first, last) / 3.0;
		for(int j = 0; j < p -> column; j++) {
			*matrix_at(p, 1, j) += *matrix_at(v, 0, j) * dist;
			*matrix_at(p, 2, j) += *matrix_at(v, 1, j) * dist;
		}
		matrix_del(v);
		addlist(ctrlp_list, p);
		return;
	}
	matrix *t = chordlength(d, first, last);
	matrix *b = matrix_bernstein(t, 3);
	matrix *p = generatecubic(d, first, last, b, v);
	int splitpoint;
	double maxerror;
	matrix *r = matrix_r_maxerror(b, p, d, first, last, maxerror, splitpoint);
	if(maxerror < error) {
		matrix_del(t); matrix_del(b); matrix_del(r); matrix_del(v);
		addlist(ctrlp_list, p);
		return;
	}
	matrix *r_deriv = r_derivative(p, t, 3);
	if(maxerror < interationerror)
		for(int i = 0; i < maxinterations; i++) {
			gauss_newton(r, r_deriv, t);
			matrix_del(r); matrix_del(r_deriv); matrix_del(p); matrix_del(b);
			b = matrix_bernstein(t, 3);
			p = generatecubic(d, first, last, b, v);
			r = matrix_r_maxerror(b, p, d, first, last, maxerror, splitpoint);
			r_deriv = r_derivative(p, t, 3);
			if(maxerror < error) {
				matrix_del(t); matrix_del(b); matrix_del(r); matrix_del(r_deriv); matrix_del(v);
				addlist(ctrlp_list, p);
				return;
			}
		}
	matrix_del(t); matrix_del(b); matrix_del(p); matrix_del(r); matrix_del(r_deriv);
	matrix *v_left = matrix_clone(v);
	matrix *v_right = matrix_clone(v);
	matrix_del(v);
	double dist = distance(d, splitpoint - 1, splitpoint + 1);
	for(int j = 0; j < d -> column; j++) {
		*matrix_at(v_left, 1, j) = (*matrix_at(d, splitpoint - 1, j) - *matrix_at(d, splitpoint + 1, j)) / dist;
		*matrix_at(v_right, 0, j) = -*matrix_at(v_left, 1, j);
	}
	fitcubic(d, first, splitpoint, v_left, ctrlp_list);
	fitcubic(d, splitpoint, last, v_right, ctrlp_list);
}

// tinh r = bp - d
matrix *matrix_r_maxerror(matrix *b, matrix *p, matrix *d, int first, int last, double &maxerror, int &splitpoint) {
	matrix *tmp = matrix_multi(b, p);
	int m = last - first + 1;
	matrix *r = matrix_create(m, d -> column);
	maxerror = -1;
	double length;
	for(int i = 0; i < m; i++) {
		length = 0;
		for(int j = 0; j < d -> column; j++) {
			*matrix_at(r, i, j) = *matrix_at(tmp, i, j) - *matrix_at(d, first + i, j);
			length += *matrix_at(r, i, j) * *matrix_at(r, i, j);
		}
		length = sqrt(length);
		if(length > maxerror) {
			maxerror = length;
			splitpoint = i + first;
		}
	}
	matrix_del(tmp);
	return r;
}

ctrlpoint_list fitcurve_cubic(matrix *d) {
	matrix *v = matrix_create(d -> column, d -> column);
	double dist1 = distance(d, 0, 1), dist2 = distance(d, d -> row - 2, d -> row - 1);
	for(int j = 0; j < d -> column; j++) {
		*matrix_at(v, 0, j) = (*matrix_at(d, 1, j) - *matrix_at(d, 0, j)) / dist1;
		*matrix_at(v, 1, j) = (*matrix_at(d, d -> row - 2, j) - *matrix_at(d, d -> row - 1, j)) / dist2;
	}
	ctrlpoint_list ctrlp_list = NULL;
	fitcubic(d, 0, d -> row - 1, v, ctrlp_list);
	return ctrlp_list;
}

void addlist(matrix_list *&list, matrix* x) {
	if(list == NULL) {
		list = (matrix_list*)malloc(sizeof(matrix_list));
		list -> mat = x;
		list -> next = NULL;
		return;
	}
	matrix_list *add = list;
	while(add -> next != NULL)
		add = add -> next;
	add -> next = (matrix_list*)malloc(sizeof(matrix_list));
	add -> next -> mat = x;
	add -> next -> next = NULL;
}

void dellist(matrix_list *&list) {
	matrix_list *tmp;
	while(list != NULL) {
		tmp = list;
		list = list -> next;
		free(tmp -> mat);
		free(tmp);
	}
	list = NULL;
}

void input(matrix *&d, int &which_method) {
	printf("[1] - Nhap du lieu tu ban phim\n[2] - Lay du lieu tu file co san\n>>");
	int kof;
	scanf("%d", &kof);
	if(kof == 1)
		getdata(d, stdin);
	else {
		int which_file;
		printf("Chon file so [1] -> [4]\n>>");
		scanf("%d", &which_file);
		FILE *f = letsopen(which_file);
		getdata(d, f);
		fclose(f);
		printf("_____________________________________________________________\n\n");
		printf("Cac diem trong file do la:\n");
		points_fprint('d', 1, d, stdout);
		printf("_____________________________________________________________\n\n");
	}
	printf("[1] - Xap xi bang 1 duong cong Bezier voi bac chi dinh\n[2] - Xap xi bang nhieu duong cong Bezier bac 3\n>>");
	scanf("%d", &which_method);
}

void getdata(matrix *&d, FILE *stream) {
	int dimension, m;
	if(stream == stdin)
		printf("Chon so chieu [2] hoac [3]:\n>>");
	fscanf(stream, "%d", &dimension);
	if(stream == stdin)
		printf("Nhap so diem:\n>>");
	fscanf(stream, "%d", &m);
	d = matrix_create(m, dimension);
	if(stream == stdin)
		printf("Nhap %d hang, moi hang la cac diem du lieu:\n", m);
	matrix_fscan(d, stream);
}

FILE *letsopen(int which_file) {
	FILE *f;
	if(which_file == 1) f = fopen("1.txt", "r");
	else if(which_file == 2) f = fopen("2.txt", "r");
	else if(which_file == 3) f = fopen("3.txt", "r");
	else f = fopen("4.txt", "r");
	return f;
}

void output(matrix_list *list, int which_method, FILE *stream) {
	fprintf(stream, "_____________________________________________________________\n\n");
	if(which_method == 1){
		fprintf(stream, "Duong cong Bezier bac %d xap xi tap diem du lieu do la:\n", list -> mat -> row - 1);
		points_fprint('p', 0, list -> mat, stream);
	}
	else {
		int count = 0;
		fprintf(stream, "Nhung duong cong Bezier bac 3 xap xi tap du lieu do la:\n");
		while(list != NULL) {
			fprintf(stream, "Duong cong thu %d:\n", ++count);
			points_fprint('p', 0, list -> mat, stream);
			list = list -> next;
		}
	}
	fprintf(stream, "_____________________________________________________________\n\n");
}