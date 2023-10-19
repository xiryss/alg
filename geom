struct point {
	ld x, y;
};
struct line {
	ld a, b, c;
	ld eval(point p) {
		return a * p.x + b * p.y + c;
	}
	void flip() {
		a = -a, b = -b, c = -c;
	}
	void norm() {
		ld d = sqrt(a * a + b * b);
		a /= d;
		b /= d;
		c /= d;
	}
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
	res.norm();
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
	res.norm();
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
ld mysqrt(ld x) {
	return (x < 0 ? 0 : sqrt(x));
}
pair<point, point> intersect(circle a, circle b) {
	if (a.r + b.r < dist(a.c, b.c)) {
		return{ {inf,inf},{inf,inf} };
	}
	if (dist(a.c, b.c) + a.r < b.r - EPS) return { {inf,inf},{inf,inf} };
	if (dist(a.c, b.c) +b.r < a.r - EPS) return { {inf,inf},{inf,inf} };
	point v = b.c - a.c;
	ld d = dist(a.c, b.c);
	ld x = (a.r * a.r - b.r * b.r + d * d) / (d * 2);
	ld ang = acos(x / a.r);
	ld cs = x / a.r;
	ld ss = mysqrt(1 - cs * cs);
	v = v * (a.r / d);
	return { a.c + spin(v, cs, ss), a.c + spin(v, cs, -ss) };
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

ld ortrsq(point a, point b) { return a * b / 2; }
ld disttoseg(point a, point s1, point s2) {
	if ((a - s1) % (s2 - s1) >= 0 && (a - s2) % (s1 - s2) >= 0) return distline(a, getln(s1, s2));
	return min(dist(a, s1), dist(a, s2));
}
bool onseg(point a, point b, point o) {
	return sgn((b - a) * (o - a)) == 0 && sgn((b - a) % (o - a)) >= 0 && sgn((a - b) % (o - b)) >= 0;
}
bool intriangle(point a, point b, point c, point o) {
	if (onseg(a, b, o) || onseg(b, c, o) || onseg(c, a, o)) return 1;
	return sgn((o - a) * (b - a)) == sgn((o - b) * (c - b)) && sgn((o - a) * (b - a)) == sgn((o - c) * (a - c));
}
vector<point> minkowski(vector<point> a, vector<point> b) {
	int pos = 0;
	for (int i = 0; i < a.size(); ++i) {
		if (a[i].y < a[pos].y || (a[i].y == a[pos].y && a[i].x < a[pos].x)) pos = i;
	}
	rotate(a.begin(), a.begin() + pos, a.end());
	pos = 0;
	for (int i = 0; i < b.size(); ++i) {
		if (b[i].y < b[pos].y || (b[i].y == b[pos].y && b[i].x < b[pos].x)) pos = i;
	}
	rotate(b.begin(), b.begin() + pos, b.end());
	int i = 0, j = 0;
	vector<point> ans;
	for (; i < a.size() || j < b.size();) {
		ans.pbc(a[i % a.size()] + b[j % b.size()]);
		ld ang = (a[(i + 1) % a.size()] - a[i % a.size()]) * (b[(j + 1) % b.size()] - b[j % b.size()]);
		if (ans.size() >= 3 && onseg(ans[ans.size() - 3], ans.back(), ans[ans.size() - 2])) {
			ans.pob;
			ans.pob;
			ans.pbc(a[i % a.size()] + b[j % b.size()]);
		}
		if (sgn(ang) > 0) ++i;
		else ++j;
	}
	if (onseg(ans[ans.size() - 2], ans[0], ans.back())) {
		ans.pob;
	}
	return ans;
}
bool inpoly(const vector<point>& a, point b) {
	int l = 1, r = a.size() - 1;
	while (r - l > 1) {
		int m = (l + r) / 2;
		if ((a[m] - a[0]) * (b - a[0]) >= 0) {
			l = m;
		}
		else r = m;
	}
	return l + 1 < a.size() && intriangle(a[0], a[l], a[l + 1], b);
}

bool bad(line a, line b, line c) {
	point x = intersect(b, c);
	assert(x.x < inf);
	return a.eval(x) > 0;
}

void tangents(point c, ld r1, ld r2, vector<line>& ans) {
	ld r = r2 - r1;
	ld z = c.x * c.x + c.y * c.y;
	ld d = z - r * r;
	if (d < -EPS)  return;
	d = sqrt(abs(d));
	line l;
	l.a = (c.x * r + c.y * d) / z;
	l.b = (c.y * r - c.x * d) / z;
	l.c = r1;
	ans.push_back(l);
}

vector<line> tangents(circle a, circle b) {
	vector<line> ans;
	for (int i = -1; i <= 1; i += 2)
		for (int j = -1; j <= 1; j += 2)
			tangents(b.c - a.c, a.r * i, b.r * j, ans);
	for (size_t i = 0; i < ans.size(); ++i)
		ans[i].c -= ans[i].a * a.c.x + ans[i].b * a.c.y;
	return ans;
}
// Do not forget about the bounding box
vector<point> hpi(vector<line> lines) {
	ld C = 1e6;
	vector<point> pts = { {-C,-C},{C,-C},{C,C},{-C,C} };
	for (auto i : lines) {
		vector<point> newpts;
		vector<point> consider;
		for (int j = 0; j < pts.size(); ++j) {
			consider.pbc(pts[j]);
			point pt = intersect(i, getln(pts[j], pts[(j + 1) % pts.size()]));
			if (onseg(pts[j], pts[(j + 1) % pts.size()], pt)) {
				consider.pbc(pt);
			}
		}
		for (auto j : consider) {
			if (i.eval(j) >= -EPS) {
				if (newpts.size() && dist(newpts.back(), j) < EPS) {
					continue;
				}
				newpts.pbc(j);
			}
		}
		if (newpts.size() && dist(newpts[0], newpts.back()) < EPS) newpts.pob;
		swap(pts, newpts);
	}
	return pts;
}
//vector<point> hpi(vector<line> lines) {
//	sort(all(lines), [](line al, line bl) -> bool {
//		point a = getvec(al);
//		point b = getvec(bl);
//		if (a.y >= 0 && b.y < 0) return 1;
//		if (a.y < 0 && b.y >= 0) return 0;
//		if (a.y == 0 && b.y == 0) return a.x > 0 && b.x < 0;
//		return (a * b) > 0;
//		});
//	for (auto i : lines) {
//		cout << i.a << ' ' << i.b << ' ' << i.c << endl;
//	}
//	vector<pair<line, int> > st;
//	for (int it = 0; it < 2; it++) {
//		for (int i = 0; i < (int)lines.size(); i++) {
//			bool flag = false;
//			while (!st.empty()) {
//				if (len(getvec(st.back().first) - getvec(lines[i])) < EPS) {
//					if (lines[i].c <= st.back().first.c) {
//						flag = true;
//						break;
//					}
//					else {
//						st.pop_back();
//					}
//				}
//				else if ((getvec(st.back().first) * getvec(lines[i])) <
//					EPS / 2) {
//					return {};
//				}
//				else if (st.size() >= 2 &&
//					bad(st[st.size() - 2].first, st[st.size() - 1].first,
//						lines[i])) {
//					st.pop_back();
//				}
//				else {
//					break;
//				}
//			}
//			if (!flag) st.push_back({ lines[i], i });
//		}
//	}
//
//	vector<int> en(lines.size(), -1);
//	vector<point> ans;
//	for (int i = 0; i < (int)st.size(); i++) {
//		if (en[st[i].second] == -1) {
//			en[st[i].second] = i;
//			continue;
//		}
//		for (int j = en[st[i].second]; j < i; j++) {
//			point I = intersect(st[j].first, st[j + 1].first);
//			ans.push_back(I);
//		}
//		break;
//	}
//	return ans;
//}
