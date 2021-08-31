struct point
{
	ld x, y;
};
struct line
{
	ld a, b, c;
};
ld det(ld a, ld b, ld c, ld d) { return a * d - b * c; }
point operator+(point a, point b) { return { a.x + b.x,a.y + b.y }; }
point operator-(point a, point b) { return { a.x - b.x,a.y - b.y }; }
point operator*(point a, ld k) { return { a.x * k, a.y * k }; }
point operator/(point a, ld k) { return { a.x / k, a.y / k }; }
line getln(point a, point b)
{
	line res;
	res.a = a.y - b.y;
	res.b = b.x - a.x;
	res.c = -(a.x * res.a + a.y * res.b);
	return res;
}
ld len(point a) { return sqrt(a.x * a.x + a.y * a.y); }
ld dist(point  a, point b) { return len(a - b); }
point intersect(line a, line b)
{
	ld zn = det(a.a, a.b, b.a, b.b);
	if (abs(zn) < EPS) return { inf,inf };
	point res;
	res.x = -det(a.c, a.b, b.c, b.b) / zn;
	res.y = -det(a.a, a.c, b.a, b.c) / zn;
	return res;
}
line perp(line a, point b)
{
	line res;
	res.a = -a.b;
	res.b = a.a;
	res.c = -(res.a * b.x + res.b * b.y);
	return res;
}
int sgn(ld x)
{
	return (x > EPS ? 1 : (x < -EPS ? -1 : 0));
}
ld val(line a, point b)
{
	return a.a * b.x + a.b * b.y + a.c;
}
