const std = @import("std");
const math = @import("std").math;
const File = @import("std").fs.File;
const io = @import("std").io;
const fs = @import("std").fs;
const rand = @import("std").rand;
const warn = @import("std").debug.warn;
const c_allocator = @import("std").heap.c_allocator;
const c_stdio = @cImport({
    // See https://github.com/zig-lang/zig/issues/515
    @cDefine("_NO_CRT_STDIO_INLINE", "1");
    @cInclude("stdio.h");
});

//float_type ft
const ft = f64;

fn Vec(comptime T: type) type {
    return struct {
        x: T,
        y: T,
        z: T,

        pub fn init(x: T, y: T, z: T) Vec(T) {
            return Vec(T){
                .x = x,
                .y = y,
                .z = z,
            };
        }

        pub fn add(self: Vec(T), other: Vec(T)) Vec(T) {
            return Vec(T).init(self.x + other.x, self.y + other.y, self.z + other.z);
        }

        pub fn sub(self: Vec(T), other: Vec(T)) Vec(T) {
            return Vec(T).init(self.x - other.x, self.y - other.y, self.z - other.z);
        }

        pub fn mul(self: Vec(T), other: Vec(T)) Vec(T) {
            return Vec(T).init(self.x * other.x, self.y * other.y, self.z * other.z);
        }

        pub fn scale(self: Vec(T), scaling: T) Vec(T) {
            return Vec(T).init(self.x * scaling, self.y * scaling, self.z * scaling);
        }

        pub fn dot(self: Vec(T), other: Vec(T)) T {
            return self.x * other.x + self.y * other.y + self.z * other.z;
        }

        pub fn cross(self: Vec(T), other: Vec(T)) Vec(T) {
            return Vec(T).init(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z, self.x * other.y - self.y * other.x);
        }

        pub fn norm(self: Vec(T)) Vec(T) {
            const factor = 1.0 / math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
            return Vec(T).init(self.x * factor, self.y * factor, self.z * factor);
        }
    };
}

const Refl_t = enum {
    DIFF,
    SPEC,
    REFR,
};

fn Ray(comptime T: type) type {
    return struct {
        o: Vec(T),
        d: Vec(T),

        pub fn init(o: Vec(T), d: Vec(T)) Ray(T) {
            return Ray(T){
                .o = o,
                .d = d,
            };
        }
    };
}

fn Sphere(comptime T: type) type {
    return struct {
        rad: T,
        p: Vec(T),
        e: Vec(T),
        c: Vec(T),
        refl: Refl_t,

        pub fn init(rad: T, p: Vec(T), e: Vec(T), c: Vec(T), refl: Refl_t) Sphere(T) {
            return Sphere(T){
                .rad = rad,
                .p = p,
                .e = e,
                .c = c,
                .refl = refl,
            };
        }

        pub fn intersect(self: Sphere(T), r: Ray(T)) T {
            const op = self.p.sub(r.o);
            const eps: T = 1e-4;
            const b = op.dot(r.d);
            var det = b * b - op.dot(op) + self.rad * self.rad;

            if (det < 0.0) {
                return 0.0;
            } else {
                det = math.sqrt(det);
            }

            const t1 = b - det;
            if (t1 > eps) {
                return t1;
            } else {
                const t2 = b + det;
                if (t2 > eps) {
                    return t2;
                } else {
                    return 0.0;
                }
            }
        }
    };
}

const spheres = [_]Sphere(ft){
    Sphere(ft).init(1.0e5, Vec(ft).init(1.0e5 + 1.0, 40.8, 81.6), Vec(ft).init(0, 0, 0), Vec(ft).init(0.75, 0.25, 0.25), Refl_t.DIFF), //Leftt
    Sphere(ft).init(1.0e5, Vec(ft).init(-1.0e5 + 99.0, 40.8, 81.6), Vec(ft).init(0, 0, 0), Vec(ft).init(0.25, 0.25, 0.75), Refl_t.DIFF), //Rght
    Sphere(ft).init(1.0e5, Vec(ft).init(50.0, 40.8, 1.0e5), Vec(ft).init(0, 0, 0), Vec(ft).init(0.75, 0.75, 0.75), Refl_t.DIFF), //Back
    Sphere(ft).init(1.0e5, Vec(ft).init(50.0, 40.8, -1.0e5 + 170.0), Vec(ft).init(0, 0, 0), Vec(ft).init(0, 0, 0), Refl_t.DIFF), //Frnt
    Sphere(ft).init(1.0e5, Vec(ft).init(50.0, 1.0e5, 81.6), Vec(ft).init(0, 0, 0), Vec(ft).init(0.75, 0.75, 0.75), Refl_t.DIFF), //Botm
    Sphere(ft).init(1.0e5, Vec(ft).init(50.0, -1.0e5 + 81.6, 81.6), Vec(ft).init(0, 0, 0), Vec(ft).init(0.75, 0.75, 0.75), Refl_t.DIFF), //Top
    Sphere(ft).init(16.5, Vec(ft).init(27.0, 16.5, 47.0), Vec(ft).init(0, 0, 0), Vec(ft).init(0.999, 0.999, 0.999), Refl_t.SPEC), //Mirr
    Sphere(ft).init(16.5, Vec(ft).init(73.0, 16.5, 78.0), Vec(ft).init(0, 0, 0), Vec(ft).init(0.999, 0.999, 0.999), Refl_t.REFR), //Glas
    Sphere(ft).init(600.0, Vec(ft).init(50.0, 681.6 - 0.27, 81.6), Vec(ft).init(12.0, 12.0, 12.0), Vec(ft).init(0, 0, 0), Refl_t.DIFF), //Lite
};

pub fn clamp(x: var) @typeOf(x) {
    if (x < 0.0) {
        return 0.0;
    } else if (x > 1.0) {
        return 1.0;
    } else {
        return x;
    }
}

pub fn toInt(x: var) i32 {
    return @floatToInt(i32, math.pow(@typeOf(x), clamp(x), 1.0 / 2.2) * 255.0 + 0.5);
}

pub fn intersect(comptime T: type, r: Ray(T), t: *T, id: *usize) bool {
    const inf = 1.0e20;
    t.* = inf;
    for (spheres) |sphere, i| {
        const d = sphere.intersect(r);
        if ((d != 0.0) and (d < t.*)) {
            t.* = d;
            id.* = i;
        }
    }

    return t.* < inf;
}

pub fn radiance(comptime T: type, r: Ray(T), depth: usize, random: *rand.Random) Vec(T) {
    var t: T = 0.0;
    var id: usize = 0;

    if (!intersect(T, r, &t, &id)) {
        return Vec(T).init(0.0, 0.0, 0.0);
    }

    const obj = spheres[id];
    const x = r.o.add(r.d.scale(t));
    const n = x.sub(obj.p).norm();
    const nl = if (n.dot(r.d) < 0) n else n.scale(-1);
    var f = obj.c;
    const p = if ((f.x > f.y) and (f.x > f.z)) f.x else if (f.y > f.z) f.y else f.z;

    if (depth > 4) {
        if ((random.float(T) < p) and depth < 10) {
            f = f.scale(1.0 / p);
        } else {
            return Vec(T).init(obj.e.x, obj.e.y, obj.e.z);
        }
    }

    if (obj.refl == Refl_t.DIFF) {
        const r1 = 2.0 * math.pi * random.float(T);
        const r2 = random.float(T);
        const r2s = math.sqrt(r2);
        const w = nl;
        const u = (if (math.fabs(w.x) > 0.1) Vec(T).init(0.0, 1.0, 0.0) else Vec(T).init(1.0, 0.0, 0.0)).cross(w).norm();
        const v = w.cross(u);
        const d = u.scale(math.cos(r1) * r2s).add(v.scale(math.sin(r1) * r2s)).add(w.scale(math.sqrt(1.0 - r2))).norm();
        return obj.e.add(f.mul(radiance(T, Ray(T).init(x, d), depth + 1, random)));
    } else if (obj.refl == Refl_t.SPEC) {
        return obj.e.add(f.mul(radiance(T, Ray(T).init(x, r.d.sub(n.scale(2.0 * n.dot(r.d)))), depth + 1, random)));
    }

    const reflRay = Ray(T).init(x, (r.d.sub(n.scale(2.0 * n.dot(r.d)))));
    const into = n.dot(nl) > 0;
    const nc: T = 1.0;
    const nt: T = 1.5;
    const nnt = if (into) nc / nt else nt / nc;
    const ddn = r.d.dot(nl);
    const cos2t = 1.0 - nnt * nnt * (1 - ddn * ddn);

    if (cos2t < 0.0) {
        return obj.e.add(f.mul(radiance(T, reflRay, depth + 1, random)));
    }

    const factor: T = if (into) 1.0 else -1.0;
    const tdir = (r.d.scale(nnt).sub(n.scale(factor * (ddn * nnt + math.sqrt(cos2t))))).norm();
    const a = nt - nc;
    const b = nt + nc;
    const c = 1.0 - (if (into) -ddn else tdir.dot(n));
    const R0 = a * a / (b * b);
    const Re = R0 + (1 - R0) * c * c * c * c * c;
    const Tr = 1.0 - Re;
    const P = 0.25 + 0.5 * Re;
    const RP = Re / P;
    const TP = Tr / (1.0 - P);

    const multi1 = if (random.float(T) < P) radiance(T, reflRay, depth + 1, random).scale(RP) else radiance(T, Ray(T).init(x, tdir), depth + 1, random).scale(TP);
    const multi2 = radiance(T, reflRay, depth + 1, random).scale(Re).add(radiance(T, Ray(T).init(x, tdir), depth + 1, random).scale(Tr));
    const multi = if (depth > 2) multi1 else multi2;

    return obj.e.add(f.mul(multi));
}

pub fn main() !void {
    var prng = rand.DefaultPrng.init(0);
    var random = prng.random;

    const w: usize = 640;
    const h: usize = 480;
    const samps = 20;

    const cam = Ray(ft).init(Vec(ft).init(50, 52, 295.6), Vec(ft).init(0, -0.042612, -1).norm());
    const cx = Vec(ft).init(@intToFloat(f64, w) * 0.5135 / @intToFloat(f64, h), 0.0, 0.0);
    const cy = (cx.cross(cam.d)).norm().scale(0.5135);
    var c = try c_allocator.alloc(Vec(ft), w * h);
    defer c_allocator.free(c);
    for (c) |*item| {
        item.* = Vec(ft).init(0.0, 0.0, 0.0);
    }
    var y: usize = 0;
    var r = Vec(ft).init(0.0, 0.0, 0.0);

    while (y < h) : (y += 1) {
        var x: usize = 0;
        while (x < w) : (x += 1) {
            const i: usize = (h - y - 1) * w + x;
            var sy: usize = 0;
            while (sy < 2) : (sy += 1) {
                var sx: usize = 0;
                while (sx < 2) : ({
                    sx += 1;
                    r = Vec(ft).init(0.0, 0.0, 0.0);
                }) {
                    var s: usize = 0;
                    while (s < samps) : (s += 1) {
                        const r1 = 2.0 * random.float(ft);
                        const dx = if (r1 < 1.0) math.sqrt(r1) - 1.0 else 1.0 - math.sqrt(2.0 - r1);
                        const r2 = 2.0 * random.float(ft);
                        const dy = if (r2 < 1.0) math.sqrt(r2) - 1.0 else 1.0 - math.sqrt(2.0 - r2);
                        const exp1 = ((@intToFloat(f64, sx) + 0.5 + dx) / 2.0 + @intToFloat(f64, x)) / @intToFloat(f64, w) - 0.5;
                        const exp2 = ((@intToFloat(f64, sy) + 0.5 + dy) / 2.0 + @intToFloat(f64, y)) / @intToFloat(f64, h) - 0.5;
                        var d = cx.scale(exp1).add(cy.scale(exp2)).add(cam.d);

                        const vec = cam.o.add(d.scale(140.0));
                        const rad = radiance(ft, Ray(ft).init(vec, d.norm()), 0, &random);
                        r = r.add(rad.scale(1.0 / @intToFloat(f64, samps)));
                    }
                    c[i] = c[i].add(Vec(ft).init(clamp(r.x), clamp(r.y), clamp(r.z))).scale(0.25);
                }
            }
        }
    }

    var file = try File.openWrite("image.ppm");
    defer file.close();
    var adapter = file.outStream();
    var buf_stream = io.BufferedOutStream(fs.File.WriteError).init(&adapter.stream);
    const stream = &buf_stream.stream;

    var i: usize = 0;
    const bits: usize = 255;
    try stream.print("P3\n{} {}\n{}\n", w, h, bits);
    while (i < (w * h)) : (i += 1) {
        try stream.print("{} {} {} ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }

    try buf_stream.flush();
}
