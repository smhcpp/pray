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
            const offset: f32 = 0.0001;
            const angle1 = corner.angle + offset;
            const angle2 = corner.angle - offset;
            const pos1 = T.Vec2f{ g.player.pos[0] + std.math.cos(angle1) * max_cast_dist, g.player.pos[1] + std.math.sin(angle1) * max_cast_dist };
            const pos2 = T.Vec2f{ g.player.pos[0] + std.math.cos(angle2) * max_cast_dist, g.player.pos[1] + std.math.sin(angle2) * max_cast_dist };
            const corner1 = Corner{ .pos = pos1, .pid = 0, .angle = angle1, .dist2 = max_dist2 };
            const corner2 = Corner{ .pos = pos2, .pid = 0, .angle = angle2, .dist2 = max_dist2 };
            _ = corner1;
            _ = corner2;
            try corners.append(g.allocator, corner);
            // try corners.append(g.allocator, corner1);
            // try corners.append(g.allocator, corner2);
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
    fn castRay(v: *Vision, start: T.Vec2f, corner: *const Corner, vision_radius_tiles: f32) ?Corner {
        const tile_size = T.WorldMap.TileSize[0];
        const vision_radius_px = vision_radius_tiles * tile_size;

        // Direction to target
        const dx = corner.pos[0] - start[0];
        const dy = corner.pos[1] - start[1];
        const ray_len = std.math.sqrt(dx * dx + dy * dy);

        if (ray_len == 0) return null;

        // Normalized direction
        const step_x = dx / ray_len;
        const step_y = dy / ray_len;

        // Current position in tile grid
        var map_x: i32 = @intFromFloat(@floor(start[0] / tile_size));
        var map_y: i32 = @intFromFloat(@floor(start[1] / tile_size));

        // Calculate step direction and initial side distances
        var step_dir_x: i32 = 1;
        var t_max_x: f32 = 0;
        var t_delta_x: f32 = 0;

        if (step_x != 0) {
            t_delta_x = @abs(tile_size / step_x);
            if (step_x > 0) {
                step_dir_x = 1;
                t_max_x = ((@as(f32, @floatFromInt(map_x)) + 1) * tile_size - start[0]) / step_x;
            } else {
                step_dir_x = -1;
                t_max_x = (start[0] - @as(f32, @floatFromInt(map_x)) * tile_size) / -step_x;
            }
        } else {
            t_delta_x = std.math.inf(f32);
            t_max_x = std.math.inf(f32);
        }

        var step_dir_y: i32 = 1;
        var t_max_y: f32 = 0;
        var t_delta_y: f32 = 0;

        if (step_y != 0) {
            t_delta_y = @abs(tile_size / step_y);
            if (step_y > 0) {
                step_dir_y = 1;
                t_max_y = ((@as(f32, @floatFromInt(map_y)) + 1) * tile_size - start[1]) / step_y;
            } else {
                step_dir_y = -1;
                t_max_y = (start[1] - @as(f32, @floatFromInt(map_y)) * tile_size) / -step_y;
            }
        } else {
            t_delta_y = std.math.inf(f32);
            t_max_y = std.math.inf(f32);
        }

        var hit_t: f32 = 0;
        var hit_horizontal = false;
        var pid: ?usize = null;

        // DDA loop
        while (hit_t < vision_radius_px and hit_t < ray_len) {
            // Check next grid cell
            if (t_max_x < t_max_y) {
                hit_t = t_max_x;
                map_x += step_dir_x;
                t_max_x += t_delta_x;
                hit_horizontal = true;
            } else {
                hit_t = t_max_y;
                map_y += step_dir_y;
                t_max_y += t_delta_y;
                hit_horizontal = false;
            }

            // Check bounds
            if (map_x < 0 or map_x >= T.WorldMap.TileNumberX or
                map_y < 0 or map_y >= T.WorldMap.TileNumberY)
            {
                pid = null;
                break;
            }

            // Check for wall
            if (v.g.wmap.tileset.get(.{ map_x, map_y })) |tile| {
                if (tile.type == .wall) {
                    pid = tile.platform_id;
                    break;
                }
            }
        }

        // If we hit the corner before hitting a wall, return the corner
        if (hit_t >= ray_len or (pid == null and hit_t >= vision_radius_px)) {
            return corner.*;
        }

        // Calculate hit position
        const hit_x = start[0] + step_x * hit_t;
        const hit_y = start[1] + step_y * hit_t;

        return Corner{
            .pos = T.Vec2f{ hit_x, hit_y },
            .pid = pid,
            .angle = corner.angle,
            .dist2 = hit_t * hit_t,
        };
    }
    fn casttRay(v: *Vision, start: T.Vec2f, corner: *const Corner, vision_radius: f32) ?Corner {
        const ray_dir = corner.pos - start;
        const tile_size = T.WorldMap.TileSize;
        const ray_len = std.math.sqrt(ray_dir[0] * ray_dir[0] + ray_dir[1] * ray_dir[1]);
        if (ray_len == 0) return null;
        const norm_ray_dir = ray_dir / T.Vec2f{ ray_len, ray_len };
        const ray_unit_distx = if (norm_ray_dir[0] != 0) @abs(1 / norm_ray_dir[0]) else std.math.inf(f32);
        const ray_unit_disty = if (norm_ray_dir[1] != 0) @abs(1 / norm_ray_dir[1]) else std.math.inf(f32);
        const ray_unit_dist = T.Vec2f{ ray_unit_distx, ray_unit_disty };
        var map_coord = @floor(start / tile_size);
        var step = std.math.sign(norm_ray_dir);
        if (step[0] == 0) step[0] = 1;
        if (step[1] == 0) step[1] = 1;
        var isHorizontal = false;
        var ray_length_bucket = T.Vec2f{ 0, 0 };
        var is_hit = false;
        var pid: ?usize = 0;
        var distance: f32 = 0;
        if (norm_ray_dir[0] < 0) {
            ray_length_bucket[0] = (start[0] / tile_size[0] - map_coord[0] + 1) * ray_unit_dist[0];
        } else {
            ray_length_bucket[0] = (map_coord[0] - start[0] / tile_size[0]) * ray_unit_dist[0];
        }
        if (norm_ray_dir[1] < 0) {
            ray_length_bucket[1] = (start[1] / tile_size[1] - map_coord[1] + 1) * ray_unit_dist[1];
        } else {
            ray_length_bucket[1] = (map_coord[1] - start[1] / tile_size[1]) * ray_unit_dist[1];
        }
        while (!is_hit) {
            if (ray_length_bucket[0] < ray_length_bucket[1]) {
                map_coord[0] += step[0];
                isHorizontal = true;
                distance = ray_length_bucket[0];
                ray_length_bucket[0] += ray_unit_dist[0];
                print("H-step: {}, distance: {}\n", .{ map_coord[0], distance });
            } else {
                map_coord[1] += step[1];
                isHorizontal = false;
                distance = ray_length_bucket[1];
                ray_length_bucket[1] += ray_unit_dist[1];
                print("V-step: {}, distance: {}\n", .{ map_coord[0], distance });
            }
            const out_of_bounds = map_coord[0] < 0 or map_coord[0] >= v.g.wmap.width or map_coord[1] < 0 or map_coord[1] >= v.g.wmap.height;
            if (out_of_bounds or distance > vision_radius) {
                is_hit = true;
                pid = null;
            }
            if (distance * tile_size[0] > ray_len) return corner.*;
            if (v.g.wmap.tileset.get(.{ @intFromFloat(map_coord[0]), @intFromFloat(map_coord[1]) })) |tile| {
                // print("hitting tile at ({}, {}), for corner pos {}\n", .{ map_coord[0], map_coord[1], corner.pos });
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
            var hit: T.Vec2f = undefined;
            if (isHorizontal) {
                hit[0] = map_coord[0] * tile_size[0]; // Exact X position of wall
                hit[1] = start[1] + norm_ray_dir[1] * dist; // Calculate Y from distance
            } else {
                hit[0] = start[0] + norm_ray_dir[0] * dist; // Calculate X from distance
                hit[1] = map_coord[1] * tile_size[1]; // Exact Y position of wall
            }
            // if (isHorizontal) {
            //     hit[0] = @round(start[0] + norm_ray_dir[0] * dist);
            //     hit[1] = start[1] + norm_ray_dir[1] * dist;
            // } else {
            //     hit[0] = start[0] + norm_ray_dir[0] * dist;
            //     hit[1] = @round(start[1] + norm_ray_dir[1] * dist);
            // }
            const epsilon = 0.01;

            // Check X
            const grid_x = @round(hit[0] / tile_size[0]) * tile_size[0];
            if (@abs(hit[0] - grid_x) < epsilon) hit[0] = grid_x;

            // Check Y
            const grid_y = @round(hit[1] / tile_size[1]) * tile_size[1];
            if (@abs(hit[1] - grid_y) < epsilon) hit[1] = grid_y;
            print("hit at ({}, {}) for corner pos {}\n", .{ hit[0], hit[1], corner.pos });
            // print("rounded hit at ({}, {})\n", .{ @round(hit[0]), @round(hit[1]) });
            // print("here is a distance: {}\n", .{dist});
            return Corner{
                .pos = hit,
                .dist2 = dist * dist,
                .angle = corner.angle,
                .pid = pid,
            };
        }
        std.process.exit(0);
        return null;
    }
    pub fn debugDrawRays(v: *Vision) void {
        const g = v.g;

        // Draw vision radius circle
        const radius_px = g.player.vision_r * T.WorldMap.TileSize[0];
        rl.drawCircleLines(@intFromFloat(g.player.pos[0]), @intFromFloat(g.player.pos[1]), radius_px, .green);

        // Draw rays to all hits
        for (v.hits.items) |hit| {
            rl.drawLineV(T.toRLVec(g.player.pos), T.toRLVec(hit.pos), .yellow);

            // Draw hit points
            rl.drawCircle(@intFromFloat(hit.pos[0]), @intFromFloat(hit.pos[1]), 3, .red);
        }
    }
    pub fn _testDrawVision(v: *Vision) void {
        if (v.hits.items.len < 3) return;
        std.mem.sort(Corner, v.hits.items, {}, comptime lessThan);
        // var write_idx: usize = 0;
        // var i: usize = 1;
        // while (i < v.hits.items.len) : (i += 1) {
        //     const current = v.hits.items[i];
        //     const prev = v.hits.items[write_idx];
        //     const angle_diff = @abs(current.angle - prev.angle);
        //     const dist_diff = @abs(current.dist2 - prev.dist2);
        //     if (angle_diff < 0.0001 and dist_diff < 1.0) {
        //         continue;
        //     }
        //     if (write_idx > 0) {
        //         const prev_prev = v.hits.items[write_idx - 1];
        //         if (current.pid == prev.pid and prev.pid == prev_prev.pid) {
        //             const p1 = prev_prev.pos;
        //             const p2 = prev.pos;
        //             const p3 = current.pos;
        //             const v1x = p2[0] - p1[0];
        //             const v1y = p2[1] - p1[1];
        //             const v2x = p3[0] - p2[0];
        //             const v2y = p3[1] - p2[1];
        //             const cross = v1x * v2y - v1y * v2x;
        //             if (@abs(cross) < 5.0) { // Tolerance of 5.0 covers float noise
        //                 v.hits.items[write_idx] = current;
        //                 continue; // Done with this point, loop again
        //             }
        //         }
        //     }
        //     write_idx += 1;
        //     v.hits.items[write_idx] = current;
        // }
        // v.hits.items.len = write_idx + 1;
        const center = v.g.player.pos;
        for (0..v.hits.items.len) |idx| {
            const a = v.hits.items[idx].pos;
            // const b = v.hits.items[(idx + 1) % v.hits.items.len].pos;
            // print("Hit Point {}: pid:{?}, dist2:{}, angle:{}\n", .{ idx, v.hits.items[idx].pid, v.hits.items[idx].dist2, v.hits.items[idx].angle });
            // rl.drawTriangle(T.toRLVec(center), T.toRLVec(b), T.toRLVec(a), .yellow);
            rl.drawLineV(T.toRLVec(center), T.toRLVec(a), .yellow);
            var buffer: [16]u8 = undefined;
            // Format 'i' into the buffer as a null-terminated string (Z)
            print("hit {} at ({}, {}) \n", .{ idx, a[0], a[1] });
            const t = std.fmt.bufPrintZ(&buffer, "{d}", .{idx}) catch " ";
            // Draw text at the hit position
            // We offset x/y slightly (-10) so the number centers on the point
            rl.drawText(t, @intFromFloat(a[0] - 10), @intFromFloat(a[1] - 20), 20, // Font Size
                .red // Color (High contrast against yellow/black)
            );
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
