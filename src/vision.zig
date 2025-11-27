const std = @import("std");
const rl = @import("raylib");
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

    pub fn updateVisionZone(v: *Vision) void {
        const g = v.g;
        const left = @max(g.player.pos[0] - g.player.vision_r, 0);
        const right = @min(g.player.pos[0] + g.player.vision_r, T.iToF32(g.screenWidth));
        const top = @max(g.player.pos[1] - g.player.vision_r, 0);
        const bottom = @min(g.player.pos[1] + g.player.vision_r, T.iToF32(g.screenHeight));
        g.wmap.platforms.items[0].pos = T.Vec2f{ left, top };
        g.wmap.platforms.items[0].size = T.Vec2f{ right - left, bottom - top };
        // std.os.exit(0);
    }

    pub fn updatePlayerVision(v: *Vision) !void {
        v.updateVisionZone();
        v.vision_step_id += 1;
        const g = v.g;
        var corners = try std.ArrayList(Corner).initCapacity(g.allocator, g.wmap.platforms.items.len);
        defer corners.deinit(g.allocator);
        const minpos = g.player.pos - T.Vec2f{ g.player.vision_r, g.player.vision_r };
        const maxpos = g.player.pos + T.Vec2f{ g.player.vision_r, g.player.vision_r };
        const minind = T.iToVec2i(minpos / T.WorldMap.GridSizeVec2f);
        const maxind = T.iToVec2i(maxpos / T.WorldMap.GridSizeVec2f);
        var j = minind[0];
        // print("pid 0: {any}\n",.{g.wmap.platforms.items[0]});
        while (j <= maxind[0]) : (j += 1) {
            var k = minind[1];
            while (k <= maxind[1]) : (k += 1) {
                if (j < 0 or k < 0 or j >= T.WorldMap.GridCellNumberX or k >= T.WorldMap.GridCellNumberY) continue;
                for (g.wmap.grid[@intCast(j)][@intCast(k)].pids.items) |pid| {
                    if (v.vision_step_id == g.wmap.platforms.items[pid].vision_step_id) continue;
                    // print("Platform ID: {}: {any}\n", .{ pid, g.wmap.platforms.items[pid] });
                    // std.process.exit(0);
                    const c = try v.getCorners(pid);
                    defer g.allocator.free(c);
                    try corners.appendSlice(g.allocator, c);
                    g.wmap.platforms.items[pid].vision_step_id = v.vision_step_id;
                }
            }
        }
        // now check pid 0 since it is not in any grid cell
        // pid 0 is the player's vision platform
        const pvis = try v.getCorners(0);
        defer g.allocator.free(pvis);
        try corners.appendSlice(g.allocator, pvis);
        try v.updateHits(&corners, minind, maxind);
    }

    fn getCorners(v: *Vision, pid: usize) ![]const Corner {
        const g = v.g;
        var corners = try std.ArrayList(Corner).initCapacity(g.allocator, 12); // Increased capacity for safety
        defer corners.deinit(g.allocator);
        const p = &g.wmap.platforms.items[pid];
        const left = @max(p.pos[0], g.player.pos[0] - g.player.vision_r, 0);
        const right = @min(p.pos[0] + p.size[0], g.player.pos[0] + g.player.vision_r, T.iToF32(g.screenWidth));
        const top = @max(p.pos[1], g.player.pos[1] - g.player.vision_r, 0);
        const bottom = @min(p.pos[1] + p.size[1], g.player.pos[1] + g.player.vision_r, T.iToF32(g.screenHeight));
        const tl = T.Vec2f{ left, top };
        const tr = T.Vec2f{ right, top };
        const bl = T.Vec2f{ left, bottom };
        const br = T.Vec2f{ right, bottom };
        const atl = std.math.atan2(tl[1] - g.player.pos[1], tl[0] - g.player.pos[0]);
        const atr = std.math.atan2(tr[1] - g.player.pos[1], tr[0] - g.player.pos[0]);
        const abl = std.math.atan2(bl[1] - g.player.pos[1], bl[0] - g.player.pos[0]);
        const abr = std.math.atan2(br[1] - g.player.pos[1], br[0] - g.player.pos[0]);
        if (left < right and top < bottom) {
            const corners_ = [_]Corner{
                Corner{ .pos = tl, .pid = pid, .angle = atl, .dist2 = T.dist2(tl, g.player.pos) },
                Corner{ .pos = tr, .pid = pid, .angle = atr, .dist2 = T.dist2(tr, g.player.pos) },
                Corner{ .pos = bl, .pid = pid, .angle = abl, .dist2 = T.dist2(bl, g.player.pos) },
                Corner{ .pos = br, .pid = pid, .angle = abr, .dist2 = T.dist2(br, g.player.pos) },
            };
            const max_cast_dist = v.g.player.vision_r * 2.0; // Cast further than vision radius
            const max_dist2 = max_cast_dist * max_cast_dist;
            for (corners_) |corner| {
                const offset: f32 = 0.0001;
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
        }
        return corners.toOwnedSlice(g.allocator);
    }

    fn updateHits(v: *Vision, corners: *std.ArrayList(Corner), minind: T.Vec2i, maxind: T.Vec2i) !void {
        if (v.vision_step_id > 10000) v.vision_step_id = 0;
        v.hits.clearRetainingCapacity();
        const g = v.g;
        for (corners.items) |corner| {
            v.vision_step_id += 1;
            var closest = Corner{ .pos = corner.pos, .pid = corner.pid, .angle = corner.angle, .dist2 = corner.dist2 };
            try v.checkVisionCollision(0, &corner, &closest);
            var i = minind[0];
            while (i <= maxind[0]) : (i += 1) {
                var j = minind[1];
                while (j <= maxind[1]) : (j += 1) {
                    if (i < 0 or j < 0 or i >= T.WorldMap.GridCellNumberX or j >= T.WorldMap.GridCellNumberY) continue;
                    for (g.wmap.grid[@intCast(i)][@intCast(j)].pids.items) |pid| {
                        if (g.wmap.platforms.items[pid].vision_step_id == v.vision_step_id) continue;
                        g.wmap.platforms.items[pid].vision_step_id = v.vision_step_id;
                        try v.checkVisionCollision(pid, &corner, &closest);
                    }
                }
            }
            try v.hits.append(v.g.allocator, closest);
        }
    }

    fn checkVisionCollision(v: *Vision, pid: usize, corner: *const Corner, closest: *Corner) !void {
        const collision = try v.getCollision(pid, corner);
        if (collision) |col| {
            const d = T.dist2(col, v.g.player.pos);
            if (closest.dist2 > d) {
                closest.pos = col;
                closest.pid = pid;
                closest.dist2 = d;
            }
        }
    }

    fn getCollision(v: *Vision, pid: usize, corner: *const Corner) !?T.Vec2f {
        const p = v.g.wmap.platforms.items[pid];
        const segtop = Segment{ .start = p.pos, .end = p.pos + T.Vec2f{ p.size[0], 0 } };
        const segleft = Segment{ .start = p.pos, .end = p.pos + T.Vec2f{ 0, p.size[1] } };
        const segright = Segment{ .start = p.pos + T.Vec2f{ p.size[0], 0 }, .end = p.pos + T.Vec2f{ p.size[0], p.size[1] } };
        const segbottom = Segment{ .start = p.pos + T.Vec2f{ 0, p.size[1] }, .end = p.pos + T.Vec2f{ p.size[0], p.size[1] } };
        const segs = [_]Segment{ segtop, segleft, segright, segbottom };
        var cols = try std.ArrayList(T.Vec2f).initCapacity(v.g.allocator, 2);
        defer cols.deinit(v.g.allocator);
        for (segs) |seg| {
            const collision = getColSegPSeg(.{ .start = v.g.player.pos, .end = corner.pos }, seg);
            if (collision) |col| {
                try cols.append(v.g.allocator, col);
            }
        }
        var closest_col: ?T.Vec2f = null;
        for (cols.items) |col| {
            if (closest_col) |closest| {
                if (T.dist2(col, v.g.player.pos) < T.dist2(closest, v.g.player.pos)) {
                    closest_col = col;
                }
            } else {
                closest_col = col;
            }
        }
        return closest_col;
    }

    pub fn drawPlayerVision(v: *Vision) void {
        std.mem.sort(Corner, v.hits.items, {}, comptime lessThan);
        var i:usize=0;
        while (i<v.hits.items.len-1):(i+=1){
            rl.drawTriangle(T.toRLVec(v.g.player.pos), T.toRLVec(v.hits.items[i+1].pos), T.toRLVec(v.hits.items[i].pos),.yellow);
        }
        rl.drawTriangle(T.toRLVec(v.g.player.pos), T.toRLVec(v.hits.items[0].pos), T.toRLVec(v.hits.items[v.hits.items.len-1].pos),.yellow);
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
