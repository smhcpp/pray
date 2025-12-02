const rl = @import("raylib");
const std = @import("std");
const Game = @import("game.zig").Game;
const math = std.math;
const print = std.debug.print;
pub const Vec2f = @Vector(2, f32);
pub const Vec2i = @Vector(2, i32);
pub fn toRLVec(vec: Vec2f) rl.Vector2 {
    return .{ .x = vec[0], .y = vec[1] };
}

pub fn dist2(a: Vec2f, b: Vec2f) f32 {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    return dx * dx + dy * dy;
}

pub fn iToVec2f(v: Vec2i) Vec2f {
    return .{ @floatFromInt(v[0]), @floatFromInt(v[1]) };
}

pub fn iToVec2i(v: Vec2f) Vec2i {
    return .{ @intFromFloat(v[0]), @intFromFloat(v[1]) };
}

pub fn iToF32(v: i32) f32 {
    return @floatFromInt(v);
}

pub fn fToI32(v: f32) i32 {
    return @intFromFloat(v);
}

/// radius is the radius of half circles on top and bottom of capsule
/// height is 4*radius and width is 2*radius
pub const Player = struct {
    pos: Vec2f,
    r: f32 = 10,
    collision_mask: u32 = 1,
    vision_mask: u32 = 1,
    // drawable: bool = true,
    texture_idle: rl.Texture2D = undefined,
    // texture_run:rl.Texture2D,
    color: rl.Color = .blue,
    vel: Vec2f = Vec2f{ 0, 0 },
    maxvel: Vec2f = Vec2f{ 400, 600 },
    speedX: f32 = 200,
    vision_r: f32 = 200,
    jump_power: f32 = 500,
    is_grounded: bool = false,

    pub fn draw(p: *Player) void {
        const baspos = toRLVec(p.pos);
        const cir1pos = baspos.subtract(toRLVec(.{ 0, p.r }));
        const cir2pos = baspos.add(toRLVec(.{ 0, p.r }));
        const recpos = baspos.subtract(toRLVec(.{ p.r, p.r }));
        rl.drawCircleV(cir1pos, p.r, p.color);
        rl.drawCircleV(cir2pos, p.r, p.color);
        rl.drawRectangleV(recpos, toRLVec(.{ p.r * 2, p.r * 2 }), p.color);
    }

    pub fn drawTexture(p: *Player) void {
        p.texture_idle.draw(fToI32(p.pos[0] - p.r), fToI32(p.pos[1] - 2 * p.r), .blue);
    }
};

/// Platforms are in the form of rectangle for now.
/// later on we can add more complex shapes like circles or polygons
/// if needed
pub const Platform = struct {
    collision_step_id: u64 = 0,
    vision_step_id: u64 = 0,
    vision_id: u32 = 1,
    collision_id: u32 = 1,
    drawable: bool = true,
    pos: Vec2f,
    size: Vec2f,
    color: rl.Color = .sky_blue,
};

pub const TileType = enum {
    empty,
    wall,
};

pub const Tile = struct {
    type: TileType = .empty,
    platform_id: usize,
};

pub fn movePlayer(g: *Game) void {
    const total_delta = g.player.vel * Vec2f{ g.dt, g.dt };
    const stepsize = Vec2f{ 2 * g.player.r - 1, 4 * g.player.r - 1 };
    const temp = total_delta / stepsize;
    const number_of_steps = @max(@floor(@abs(temp[0])), @floor(@abs(temp[1]))) + 1;
    var bucket = Vec2f{ 0, 0 };
    var poschange = Vec2f{ 0, 0 };
    var velocity_step = total_delta / Vec2f{ number_of_steps, number_of_steps };
    var i: f32 = 0;
    var retvel = Vec2f{ 0, 0 };
    var xcollision = false;
    var ycollision = false;
    while ((!xcollision or !ycollision) and i < number_of_steps) : (i += 1) {
        g.collision_step_id += 1;
        bucket += velocity_step;
        const pos = g.player.pos + bucket;
        const minpos = pos - Vec2f{ g.player.r, 2 * g.player.r };
        const maxpos = pos + Vec2f{ g.player.r, 2 * g.player.r };
        const minind = iToVec2i(minpos / WorldMap.GridSizeVec2f);
        const maxind = iToVec2i(maxpos / WorldMap.GridSizeVec2f);
        var j = minind[0];
        var k = minind[1];
        while (j <= maxind[0]) : (j += 1) {
            while (k <= maxind[1]) : (k += 1) {
                if (j < 0 or k < 0 or j >= WorldMap.GridCellNumberX or k >= WorldMap.GridCellNumberY) continue;
                for (g.wmap.grid[@intCast(j)][@intCast(k)].pids.items) |pid| {
                    if (g.collision_step_id == g.wmap.platforms.items[pid].collision_step_id) continue;
                    const col_id = g.wmap.platforms.items[pid].collision_id;
                    if (g.player.collision_mask | col_id != col_id and col_id == 0) continue;
                    retvel = checkPlayerCollision(pos, g.player.r, g.player.vel, g.wmap.platforms.items[pid], &poschange);
                    if (g.player.vel[0] != retvel[0]) {
                        velocity_step[0] = 0;
                        xcollision = true;
                    }
                    if (g.player.vel[1] != retvel[1]) {
                        ycollision = true;
                        velocity_step[1] = 0;
                    }
                    g.wmap.platforms.items[pid].collision_step_id = g.collision_step_id;
                }
            }
        }
    }
    if (ycollision) retvel[1] = 0;
    if (xcollision) retvel[0] = 0;
    g.player.pos += bucket + poschange;
    if (xcollision or ycollision) g.player.vel = retvel;
}

pub fn checkPlayerCollision(pos: Vec2f, capr: f32, vel: Vec2f, plat: Platform, poschange: *Vec2f) Vec2f {
    var retvel = vel;
    const captop = pos[1] - capr;
    const capbot = pos[1] + capr;
    const closest_plat_point_x = math.clamp(pos[0], plat.pos[0], plat.pos[0] + plat.size[0]);
    const closest_plat_point_y = math.clamp(pos[1], plat.pos[1], plat.pos[1] + plat.size[1]);
    const closest_cap_point_y = math.clamp(closest_plat_point_y, captop, capbot);
    const closest_cap_point_x = pos[0];
    const dx = closest_cap_point_x - closest_plat_point_x;
    const dy = closest_cap_point_y - closest_plat_point_y;
    const dis2 = dx * dx + dy * dy;
    if (dis2 < capr * capr and dis2 > 0) {
        const dist = math.sqrt(dis2);
        const overlap = capr - dist;
        const dx_norm = dx / dist;
        const dy_norm = dy / dist;
        poschange[0] += dx_norm * overlap;
        poschange[1] += dy_norm * overlap;
        if (dy_norm < -0.7) {
            if (vel[1] > 0) retvel[1] = 0;
        }
        if (dy_norm > 0.7 and retvel[1] < 0) {
            retvel[1] = 0;
        }
        if (@abs(dx_norm) > 0.7) {
            retvel[0] = 0;
        }
    }
    return retvel;
}

pub const WorldMap = struct {
    pub const TileSize: f32 = 8;
    pub const TileNumberX: f32 = 128;
    pub const TileNumberY: f32 = 100;
    pub const TileSizeVec2f = Vec2f{ TileSize, TileSize };
    tileset: std.AutoHashMap(Vec2i, Tile),
    platforms: std.ArrayList(Platform),

    pub fn init(allocator: std.mem.Allocator) !*WorldMap {
        const map = try allocator.create(WorldMap);
        map.* = .{
            .platforms = try std.ArrayList(Platform).initCapacity(allocator, 10),
            .tileset = std.AutoHashMap(Vec2i, Tile).init(allocator),
        };
        try map.setup(allocator);
        return map;
    }

    pub fn setup(map: *WorldMap, allocator: std.mem.Allocator) !void {
        try map.platforms.append(allocator, Platform{
            .pos = .{ 0, 0 },
            .size = .{ 1, 1 },
            .collision_id = 0,
            .drawable = false,
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 13, 6 },
            .size = .{ 12, 24 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 24, 24 },
            .size = .{ 36, 18 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 85, 28 },
            .size = .{ 24, 18 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 60, 54 },
            .size = .{ 24, 12 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 0, TileNumberY - 1 },
            .size = .{  TileNumberX , 1 },
        });

        for (map.platforms.items, 0..) |platform, pid| {
            if (platform.collision_id == 0) continue;
            const imin: usize = @intCast(platform.pos[0]);
            const imax: usize = @intCast(platform.pos[0] + platform.size[0]);
            const jmin: usize = @intCast(platform.pos[1]);
            const jmax: usize = @intCast(platform.pos[1] + platform.size[1]);
            var i = imin;
            while (i <= imax) : (i += 1) {
                var j = jmin;
                while (j <= jmax) : (j += 1) {
                    map.tileset.put(.{ @intCast(i), @intCast(j) }, .{ .type = .block, .platform_id = pid });
                }
            }
        }
    }

    pub fn deinit(map: *WorldMap, allocator: std.mem.Allocator) void {
        map.platforms.deinit(allocator);
        map.tileset.deinit();
        allocator.destroy(map);
    }
};
