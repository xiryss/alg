struct point {
	ld x, y;
};
struct line {
	ld a, b, c;
};
ld det(ld a, ld b, ld c, ld d) { return a * d - b * c; }
point operator+(point a, point b) { return { a.x + b.x,a.y + b.y }; }
point operator-(point a, point b) { return { a.x - b.x,a.y - b.y }; }
point operator*(point a, ld k) { return { a.x * k, a.y * k }; }
point operator/(point a, ld k) { return { a.x / k, a.y / k }; }
ld operator*(point a, point b) { return a.x * b.y - a.y * b.x; }
ld operator%(point a, point b) { return a.x * b.x + a.y * b.y; }
line getln(point a, point b) {
	line res;
	res.a = a.y - b.y;
	res.b = b.x - a.x;
	res.c = -(a.x * res.a + a.y * res.b);
	return res;
}
ld len(point a) { return sqrt(a.x * a.x + a.y * a.y); }
ld dist(point  a, point b) { return len(a - b); }
point intersect(line a, line b) {
	ld zn = det(a.a, a.b, b.a, b.b);
	if (abs(zn) < EPS) return { inf,inf };
	point res;
	res.x = -det(a.c, a.b, b.c, b.b) / zn;
	res.y = -det(a.a, a.c, b.a, b.c) / zn;
	return res;
}
line perp(line a, point b) {
	line res;
	res.a = -a.b;
	res.b = a.a;
	res.c = -(res.a * b.x + res.b * b.y);
	return res;
}
int sgn(ld x) {
	return (x > EPS ? 1 : (x < -EPS ? -1 : 0));
}
ld val(line a, point b) {
	return a.a * b.x + a.b * b.y + a.c;
}
point spin(point a, ld cs, ld ss) {
	return { a.x * cs - a.y * ss, a.x * ss + a.y * cs };
}
struct circle {
	point c;
	ld r;
};
point getvec(line a) {
	return { -a.b, a.a };
}
pair<point, point> intersect(circle a, line b) {
	line pp = perp(b, a.c);
	point xx = intersect(pp, b);
	ld ds = dist(a.c, xx);
	ld td = sqrt(a.r * a.r - ds * ds);
	point vv = getvec(b); vv = vv / len(vv) * td;
	return { xx + vv, xx - vv };
}
ld distline(point a, line b) {
	return abs(val(b, a)) / sqrt(b.a * b.a + b.b * b.b);
}
ld getangle(point a, point b) {
	return atan2(a * b, a % b);
}
ld trsq(point a, point b) {
	return abs(a * b) / 2;
}
ld disttoseg(point a, point s1, point s2) {
	if ((a - s1) % (s2 - s1) >= 0 && (a - s2) % (s1 - s2) >= 0) return distline(a, getln(s1, s2));
	return min(dist(a, s1), dist(a, s2));
}
