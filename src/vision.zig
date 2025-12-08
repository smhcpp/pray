const std = @import("std");
const rl = @import("raylib");
const gl = rl.gl;
const T = @import("types.zig");
const print = std.debug.print;
const Game = @import("game.zig").Game;
pub const Vision = struct {
    vision_step_id: u64 = 0,
    hits: std.ArrayList(Corner),
    g: *Game,

    pub fn init(g: *Game) !*Vision {
        const v = try g.allocator.create(Vision);
        v.* = Vision{
            .g = g,
            .hits = try std.ArrayList(Corner).initCapacity(g.allocator, g.wmap.platforms.items.len),
        };
        return v;
    }

    pub fn deinit(v: *Vision) void {
        v.hits.deinit(v.g.allocator);
        v.g.allocator.destroy(v);
    }

    fn insertVisionCorners(v: *Vision, corners: *std.ArrayList(Corner)) !void {
        const g = v.g;
        const tile_size = T.WorldMap.TileSize[0];
        const left = @max(g.player.pos[0] - g.player.vision_r * tile_size, 0);
        const top = @max(g.player.pos[1] - g.player.vision_r * tile_size, 0);
        const right = @min(g.player.pos[0] + g.player.vision_r * tile_size, T.iToF32(g.screenWidth));
        const bottom = @min(g.player.pos[1] + g.player.vision_r * tile_size, T.iToF32(g.screenHeight));
        if (left < right and top < bottom) {
            const c = try v.getCorners(.{ top, right, bottom, left }, null);
            defer g.allocator.free(c);
            try corners.appendSlice(g.allocator, c);
        }
    }

    pub fn updatePlayerVision(v: *Vision) !void {
        const tile_size = T.WorldMap.TileSize[0];
        const g = v.g;
        var corners = try std.ArrayList(Corner).initCapacity(g.allocator, g.wmap.platforms.items.len);
        defer corners.deinit(g.allocator);
        try v.insertVisionCorners(&corners);
        for (g.wmap.platforms.items, 0..) |p, pid| {
            const left = @max(p.pos[0] * tile_size, g.player.pos[0] - g.player.vision_r * tile_size, 0);
            const top = @max(p.pos[1] * tile_size, g.player.pos[1] - g.player.vision_r * tile_size, 0);
            const right = @min(p.pos[0] * tile_size + p.size[0] * tile_size, g.player.pos[0] + g.player.vision_r * tile_size, T.iToF32(g.screenWidth));
            const bottom = @min(p.pos[1] * tile_size + p.size[1] * tile_size, g.player.pos[1] + g.player.vision_r * tile_size, T.iToF32(g.screenHeight));
            if (left < right and top < bottom) {
                // print("Platform ID: {}: {any}\n", .{ pid, g.wmap.platforms.items[pid] });
                // std.process.exit(0);
                const c = try v.getCorners(.{ top, right, bottom, left }, pid);
                defer g.allocator.free(c);
                try corners.appendSlice(g.allocator, c);
            }
        }
        try v.updateHits(&corners);
    }

    fn getCorners(v: *Vision, sides: [4]f32, pid: ?usize) ![]const Corner {
        const g = v.g;
        const tile_size = T.WorldMap.TileSize[0];
        var corners = try std.ArrayList(Corner).initCapacity(g.allocator, 12); // Increased capacity for safety
        defer corners.deinit(g.allocator);
        const tl = T.Vec2f{ sides[3], sides[0] };
        const tr = T.Vec2f{ sides[1], sides[0] };
        const bl = T.Vec2f{ sides[3], sides[2] };
        const br = T.Vec2f{ sides[1], sides[2] };
        const atl = std.math.atan2(tl[1] - g.player.pos[1], tl[0] - g.player.pos[0]);
        const atr = std.math.atan2(tr[1] - g.player.pos[1], tr[0] - g.player.pos[0]);
        const abl = std.math.atan2(bl[1] - g.player.pos[1], bl[0] - g.player.pos[0]);
        const abr = std.math.atan2(br[1] - g.player.pos[1], br[0] - g.player.pos[0]);
        const corners_ = [_]Corner{
            Corner{ .pos = tl, .pid = pid, .angle = atl, .dist2 = T.dist2(tl, g.player.pos) },
            Corner{ .pos = tr, .pid = pid, .angle = atr, .dist2 = T.dist2(tr, g.player.pos) },
            Corner{ .pos = bl, .pid = pid, .angle = abl, .dist2 = T.dist2(bl, g.player.pos) },
            Corner{ .pos = br, .pid = pid, .angle = abr, .dist2 = T.dist2(br, g.player.pos) },
        };
        const max_cast_dist = v.g.player.vision_r * tile_size; // Cast further than vision radius
        const max_dist2 = max_cast_dist * max_cast_dist * 2;
        for (corners_) |corner| {
            const offset: f32 = 0.01;
            const angle1 = corner.angle + offset;
            const angle2 = corner.angle - offset;
            const pos1 = T.Vec2f{ g.player.pos[0] + std.math.cos(angle1) * max_cast_dist, g.player.pos[1] + std.math.sin(angle1) * max_cast_dist };
            const pos2 = T.Vec2f{ g.player.pos[0] + std.math.cos(angle2) * max_cast_dist, g.player.pos[1] + std.math.sin(angle2) * max_cast_dist };
            const corner1 = Corner{ .pos = pos1, .pid = 0, .angle = angle1, .dist2 = max_dist2 };
            const corner2 = Corner{ .pos = pos2, .pid = 0, .angle = angle2, .dist2 = max_dist2 };
            try corners.append(g.allocator, corner);
            try corners.append(g.allocator, corner1);
            try corners.append(g.allocator, corner2);
        }
        return corners.toOwnedSlice(g.allocator);
    }

    fn updateHits(v: *Vision, corners: *std.ArrayList(Corner)) !void {
        if (v.vision_step_id > 10000) v.vision_step_id = 0;
        v.hits.clearRetainingCapacity();
        const g = v.g;
        print("--------------------------------\n", .{});
        for (corners.items) |corner| {
            var closest = Corner{ .pos = corner.pos, .pid = corner.pid, .angle = corner.angle, .dist2 = corner.dist2 };
            v.checkVisionCollision(&corner, &closest);
            // try v.checkPenetratedVisionCollision( &corner, &closest);
            try v.hits.append(g.allocator, closest);
        }
    }

    fn checkVisionCollision(v: *Vision, corner: *const Corner, closest: *Corner) void {
        const collision = v.castRay(v.g.player.pos, corner, v.g.player.vision_r);
        if (collision) |col| {
            const d = T.dist2(col.pos, v.g.player.pos);
            if (closest.dist2 > d) {
                closest.pos = col.pos;
                closest.pid = col.pid;
                closest.dist2 = col.dist2;
            }
        }
    }

    fn castRay(v: *Vision, start: T.Vec2f, corner: *const Corner, vision_radius: f32) ?Corner {
        const ray_dir = corner.pos - start;
        const tile_size = T.WorldMap.TileSize;
        const ray_len = std.math.sqrt(ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1]);
        if (ray_len == 0) return null;
        const norm_ray_dir = ray_dir / T.Vec2f{ ray_len, ray_len };
        const ray_unit_distx = if (norm_ray_dir[0] != 0) @abs(1 / norm_ray_dir[0]) else std.math.inf(f32);
        const ray_unit_disty = if (norm_ray_dir[1] != 0) @abs(1 / norm_ray_dir[1]) else std.math.inf(f32);
        const ray_unit_dist = T.Vec2f{ ray_unit_distx, ray_unit_disty };
        var map_coord = @floor(start / tile_size);
        const start_coord = map_coord;
        var step = T.Vec2f{ 0, 0 };
        var isHorizontal = false;
        var ray_length_bucket = T.Vec2f{ 0, 0 };
        var is_hit = false;
        var pid: ?usize = 0;
        if (norm_ray_dir[0] < 0) {
            step[0] = -1;
            ray_length_bucket[0] = (start[0] / tile_size[0] - map_coord[0] + 1) * ray_unit_dist[0];
        } else {
            step[0] = 1;
            ray_length_bucket[0] = (map_coord[0] + 1 - start[0] / tile_size[0]) * ray_unit_dist[0];
        }
        if (norm_ray_dir[1] < 0) {
            step[1] = -1;
            ray_length_bucket[1] = (start[1] / tile_size[1] - map_coord[1] + 1) * ray_unit_dist[1];
        } else {
            step[1] = 1;
            ray_length_bucket[1] = (map_coord[1] + 1 - start[1] / tile_size[1]) * ray_unit_dist[1];
        }
        while (!is_hit) {
            if (ray_length_bucket[0] < ray_length_bucket[1]) {
                map_coord[0] += step[0];
                isHorizontal = true;
                ray_length_bucket[0] += ray_unit_dist[0];
            } else {
                map_coord[1] += step[1];
                isHorizontal = false;
                ray_length_bucket[1] += ray_unit_dist[1];
            }
            const distance = @abs(map_coord - start_coord);
            const out_of_bounds = map_coord[0] < 0 or map_coord[0] >= v.g.wmap.width or map_coord[1] < 0 or map_coord[1] >= v.g.wmap.height;
            if (out_of_bounds or distance[0] > vision_radius or distance[1] > vision_radius) {
                is_hit = true;
                pid = null;
            }
            if (v.g.wmap.tileset.get(.{ @intFromFloat(map_coord[0]), @intFromFloat(map_coord[1]) })) |tile| {
                if (tile.type == .wall) {
                    is_hit = true;
                    pid = tile.platform_id;
                }
            }
        }
        if (is_hit) {
            const realdist = ray_length_bucket - ray_unit_dist;
            const dist = if (isHorizontal) realdist[0] * tile_size[0] else realdist[1] * tile_size[1];
            if (dist * dist > corner.dist2) return corner.*;
            const hit = start + norm_ray_dir * T.Vec2f{ dist, dist };
            // print("here is a distance: {}\n", .{dist});
            return Corner{
                .pos = hit,
                .dist2 = dist * dist,
                .angle = corner.angle,
                .pid = pid,
            };
        }
        return null;
    }

    pub fn _testDrawVision(v: *Vision) void {
        if (v.hits.items.len < 3) return;
        std.mem.sort(Corner, v.hits.items, {}, comptime lessThan);
        var write_idx: usize = 0;
        var i: usize = 1;
        while (i < v.hits.items.len) : (i += 1) {
            const current = v.hits.items[i];
            const prev = v.hits.items[write_idx];
            const angle_diff = @abs(current.angle - prev.angle);
            const dist_diff = @abs(current.dist2 - prev.dist2);
            if (angle_diff < 0.0001 and dist_diff < 1.0) {
                continue;
            }
            if (write_idx > 0) {
                const prev_prev = v.hits.items[write_idx - 1];
                if (current.pid == prev.pid and prev.pid == prev_prev.pid) {
                    const p1 = prev_prev.pos;
                    const p2 = prev.pos;
                    const p3 = current.pos;
                    const v1x = p2[0] - p1[0];
                    const v1y = p2[1] - p1[1];
                    const v2x = p3[0] - p2[0];
                    const v2y = p3[1] - p2[1];
                    const cross = v1x * v2y - v1y * v2x;
                    if (@abs(cross) < 5.0) { // Tolerance of 5.0 covers float noise
                        v.hits.items[write_idx] = current;
                        continue; // Done with this point, loop again
                    }
                }
            }
            write_idx += 1;
            v.hits.items[write_idx] = current;
        }
        v.hits.items.len = write_idx + 1;
        const center = v.g.player.pos;
        for (0..v.hits.items.len) |idx| {
            const a = v.hits.items[idx].pos;
            const b = v.hits.items[(idx + 1) % v.hits.items.len].pos;
            print("Hit Point {}: pid:{?}, dist2:{}, angle:{}\n", .{ idx, v.hits.items[idx].pid, v.hits.items[idx].dist2, v.hits.items[idx].angle });
            rl.drawTriangle(T.toRLVec(center), T.toRLVec(b), T.toRLVec(a), .yellow);
        }
    }

    pub fn drawPlayerVision(v: *Vision) void {
        if (v.hits.items.len < 3) return;
        std.mem.sort(Corner, v.hits.items, {}, comptime lessThan);
        const center = v.g.player.pos;
        gl.rlBegin(gl.rl_triangles);
        for (0..v.hits.items.len) |i| {
            print("Hit Point {}: dist2:{} and angle:{}\n", .{ i, v.hits.items[i].dist2, v.hits.items[i].angle });
            const a = v.hits.items[i].pos;
            const b = v.hits.items[(i + 1) % v.hits.items.len].pos;
            gl.rlColor4ub(255, 255, 255, 255);
            gl.rlVertex2f(center[0], center[1]);

            gl.rlColor4ub(255, 255, 255, 1);
            gl.rlVertex2f(b[0], b[1]);
            gl.rlVertex2f(a[0], a[1]);
        }

        gl.rlEnd();
    }
};

fn lessThan(context: void, a: Corner, b: Corner) bool {
    _ = context;
    return a.angle < b.angle;
}
// Returns intersection point if 'seg' intersects 'pseg'.
// Assumes 'pseg' is strictly vertical or horizontal.
pub fn getColSegPSeg(seg: Segment, pseg: Segment) ?T.Vec2f {
    const is_horizontal = pseg.start[1] == pseg.end[1];
    if (is_horizontal) {
        const wall_y = pseg.start[1];
        if ((seg.start[1] > wall_y and seg.end[1] > wall_y) or
            (seg.start[1] < wall_y and seg.end[1] < wall_y))
        {
            return null;
        }
        const dy = seg.end[1] - seg.start[1];
        if (@abs(dy) < 1e-3) return null;
        const t = (wall_y - seg.start[1]) / dy;
        if (t < 0 or t > 1) return null;
        const intersect_x = seg.start[0] + t * (seg.end[0] - seg.start[0]);
        const min_x = @min(pseg.start[0], pseg.end[0]);
        const max_x = @max(pseg.start[0], pseg.end[0]);
        if (intersect_x >= min_x and intersect_x <= max_x) {
            return T.Vec2f{ intersect_x, wall_y };
        }
    } else {
        const wall_x = pseg.start[0];
        if ((seg.start[0] > wall_x and seg.end[0] > wall_x) or
            (seg.start[0] < wall_x and seg.end[0] < wall_x))
        {
            return null;
        }
        const dx = seg.end[0] - seg.start[0];
        if (@abs(dx) < 1e-6) return null;
        const t = (wall_x - seg.start[0]) / dx;
        if (t < 0 or t > 1) return null;
        const intersect_y = seg.start[1] + t * (seg.end[1] - seg.start[1]);
        const min_y = @min(pseg.start[1], pseg.end[1]);
        const max_y = @max(pseg.start[1], pseg.end[1]);
        if (intersect_y >= min_y and intersect_y <= max_y) {
            return T.Vec2f{ wall_x, intersect_y };
        }
    }
    return null;
}

pub const Segment = struct {
    start: T.Vec2f,
    end: T.Vec2f,
};

pub const Corner = struct {
    pos: T.Vec2f,
    pid: ?usize,
    angle: f32,
    dist2: f32,
};
