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
    color: rl.Color = .blue,
    vel: Vec2f = Vec2f{ 0, 0 },
    maxvel: Vec2f = Vec2f{ 400, 600 },
    speedX: f32 = 250,
    vision_r: f32 = 200,
    jump_power: f32 = 500,
    is_grounded: bool = false,

    pub fn draw(g: *Player) void {
        const baspos = toRLVec(g.pos);
        const cir1pos = baspos.subtract(toRLVec(.{ 0, g.r }));
        const cir2pos = baspos.add(toRLVec(.{ 0, g.r }));
        const recpos = baspos.subtract(toRLVec(.{ g.r, g.r }));
        rl.drawCircleV(cir1pos, g.r, g.color);
        rl.drawCircleV(cir2pos, g.r, g.color);
        rl.drawRectangleV(recpos, toRLVec(.{ g.r * 2, g.r * 2 }), g.color);
    }
};

// pub const VisionMask = enum(u32) {
// };

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
    color: rl.Color = .gray,
};

pub const GridCell = struct {
    pids: std.ArrayList(usize),
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
    // find closest point on platform to the capsule's center
    const closest_plat_point_x = math.clamp(pos[0], plat.pos[0], plat.pos[0] + plat.size[0]);
    const closest_plat_point_y = math.clamp(pos[1], plat.pos[1], plat.pos[1] + plat.size[1]);

    // find closest point on capsule central vertical segment to the closest point of the platform
    // that we found above
    const closest_cap_point_y = math.clamp(closest_plat_point_y, captop, capbot);
    const closest_cap_point_x = pos[0];

    // find their distance^2:
    const dx = closest_cap_point_x - closest_plat_point_x;
    const dy = closest_cap_point_y - closest_plat_point_y;
    const dis2 = dx * dx + dy * dy;

    // if there is collision
    if (dis2 < capr * capr and dis2 > 0) {
        const dist = math.sqrt(dis2);
        const overlap = capr - dist;
        const dx_norm = dx / dist;
        const dy_norm = dy / dist;

        poschange[0] += dx_norm * overlap;
        poschange[1] += dy_norm * overlap;

        // if player is moving down and colliding with platform
        if (dy_norm < -0.7) {
            if (vel[1] > 0) retvel[1] = 0;
        }
        // if player is hitting the platform from below
        if (dy_norm > 0.7 and retvel[1] < 0) {
            retvel[1] = 0;
        }
        // if player is hitting the platform from the side
        if (@abs(dx_norm) > 0.7) {
            retvel[0] = 0;
        }
    }
    return retvel;
}

pub const WorldMap = struct {
    pub const GridSize: f32 = 128;
    pub const GridCellNumberX: f32 = 10;
    pub const GridCellNumberY: f32 = 6;
    pub const GridSizeVec2f = Vec2f{ GridSize, GridSize };
    grid: [GridCellNumberX][GridCellNumberY]GridCell,
    platforms: std.ArrayList(Platform),

    pub fn init(allocator: std.mem.Allocator) !*WorldMap {
        const map = try allocator.create(WorldMap);
        map.* = .{
            .grid = undefined,
            .platforms = try std.ArrayList(Platform).initCapacity(allocator, 10),
        };
        try map.setup(allocator);
        return map;
    }

    pub fn setup(map: *WorldMap, allocator: std.mem.Allocator) !void {
        for (0..GridCellNumberX) |x| {
            for (0..GridCellNumberY) |y| {
                map.grid[x][y] = .{ .pids = try std.ArrayList(usize).initCapacity(allocator, 10) };
            }
        }

        try map.platforms.append(allocator, Platform{
            .pos = .{ 0, 0 },
            .size = .{1,1},
            .collision_id = 0,
            .drawable = false,
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 100, 50 },
            .size = .{ 100, 200 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 200, 200 },
            .size = .{ 300, 150 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 700, 250 },
            .size = .{ 200, 150 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 500, 450 },
            .size = .{ 200, 100 },
        });

        try map.platforms.append(allocator, Platform{
            .pos = .{ 0, GridSize * GridCellNumberY - 11 },
            .size = .{ GridSize * GridCellNumberX - 1, 10 },
        });

        for (map.platforms.items, 0..) |platform, pid| {
            if (pid == 0) continue;
            const imin: usize = @intFromFloat(platform.pos[0] / iToF32(WorldMap.GridSize));
            const jmin: usize = @intFromFloat(platform.pos[1] / iToF32(WorldMap.GridSize));
            const imax: usize = @intFromFloat((platform.pos[0] + platform.size[0]) / iToF32(WorldMap.GridSize));
            const jmax: usize = @intFromFloat((platform.pos[1] + platform.size[1]) / iToF32(WorldMap.GridSize));
            var i = imin;
            while (i <= imax) : (i += 1) {
                var j = jmin;
                while (j <= jmax) : (j += 1) {
                    try map.grid[i][j].pids.append(allocator, pid);
                }
            }
        }
    }

    pub fn deinit(map: *WorldMap, allocator: std.mem.Allocator) void {
        for (0..GridCellNumberX) |x| {
            for (0..GridCellNumberY) |y| {
                map.grid[x][y].pids.deinit(allocator);
            }
        }
        map.platforms.deinit(allocator);
        allocator.destroy(map);
    }
};
