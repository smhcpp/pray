const rl = @import("raylib");
const std = @import("std");
const T = @import("types.zig");
const Vision = @import("vision.zig").Vision;
const math = std.math;
const print = std.debug.print;
pub const Game = struct {
    collision_step_id: u64 = 0,
    vision: *Vision = undefined,
    pause: bool = false,
    fps: i32 = 60,
    gravity: f32 = 500,
    friction: f32 = 10,
    dt: f32 = undefined,
    allocator: std.mem.Allocator,
    player: T.Player = undefined,
    wmap: *T.WorldMap = undefined,
    screenWidth: i32 = T.fToI32(T.WorldMap.GridSize * T.WorldMap.GridCellNumberX),
    screenHeight: i32 = T.fToI32(T.WorldMap.GridSize * T.WorldMap.GridCellNumberY),
    inputs: struct {
        right: bool,
        left: bool,
        up: bool,
        down: bool,
        jump: bool,
        dash: bool,
        attack: bool,
        escape: bool,
    } = .{
        .right = false,
        .left = false,
        .up = false,
        .down = false,
        .jump = false,
        .dash = false,
        .attack = false,
        .escape = false,
    },

    pub fn init(allocator: std.mem.Allocator) !*Game {
        const game = try allocator.create(Game);
        game.* = .{
            .allocator = allocator,
        };
        try game.setup();
        game.dt = 1.0 / T.iToF32(game.fps);
        return game;
    }

    fn _testInput(g: *Game) void {
        const delta: f32 = 3;
        if (g.inputs.right) {
            g.player.pos[0] += delta;
        } else if (g.inputs.left) {
            g.player.pos[0] -= delta;
        } else if (g.inputs.up) {
            g.player.pos[1] -= delta;
        } else if (g.inputs.down) {
            g.player.pos[1] += delta;
        }
    }

    fn updateInputs(g: *Game) void {
        g.inputs.right = rl.isKeyDown(rl.KeyboardKey.right);
        g.inputs.left = rl.isKeyDown(rl.KeyboardKey.left);
        g.inputs.up = rl.isKeyDown(rl.KeyboardKey.up);
        g.inputs.down = rl.isKeyDown(rl.KeyboardKey.down);
        g.inputs.jump = rl.isKeyDown(rl.KeyboardKey.space);
        g.inputs.dash = rl.isKeyDown(rl.KeyboardKey.e);
        g.inputs.attack = rl.isKeyDown(rl.KeyboardKey.r);
        g.inputs.escape = rl.isKeyPressed(rl.KeyboardKey.escape);
        if (g.inputs.escape) {
            g.pause = !g.pause;
            print("Pause: {}\n", .{g.pause});
        }
    }

    fn setup(g: *Game) !void {
        g.wmap = try T.WorldMap.init(g.allocator);
        g.vision = try Vision.init(g);
        g.dt = 1.0 / T.iToF32(g.fps);
        g.player = T.Player{
            .pos = .{ T.iToF32(@divTrunc(g.screenWidth, 2)), T.iToF32(@divTrunc(g.screenHeight, 2)) },
        };
    }

    fn draw(g: *Game) void {
        if (g.pause) g.drawGridLines();
        for (g.wmap.platforms.items) |p| {
            if (!p.drawable) continue;
            const left = @max(p.pos[0], g.player.pos[0] - g.player.vision_r, 0);
            const right = @min(p.pos[0] + p.size[0], g.player.pos[0] + g.player.vision_r, T.iToF32(g.screenWidth));
            const top = @max(p.pos[1], g.player.pos[1] - g.player.vision_r, 0);
            const bottom = @min(p.pos[1] + p.size[1], g.player.pos[1] + g.player.vision_r, T.iToF32(g.screenHeight));
            if (left < right and top < bottom)
                rl.drawRectangleV(.{ .x = left, .y = top }, .{ .x = right - left, .y = bottom - top }, p.color);
        }
        g.vision.drawPlayerVision();
        g.player.draw();
    }

    fn drawGridLines(g: *Game) void {
        var i: i32 = 0;
        while (i < T.WorldMap.GridCellNumberY) : (i += 1) {
            rl.drawLineV(rl.Vector2{ .x = 0, .y = T.iToF32(i) * T.WorldMap.GridSize }, rl.Vector2{ .x = T.iToF32(g.screenWidth), .y = T.iToF32(i) * T.WorldMap.GridSize }, .yellow);
        }
        i = 0;
        while (i < T.WorldMap.GridCellNumberX) : (i += 1) {
            rl.drawLineV(rl.Vector2{ .x = T.iToF32(i) * T.WorldMap.GridSize, .y = 0 }, rl.Vector2{ .x = T.iToF32(i) * T.WorldMap.GridSize, .y = T.iToF32(g.screenHeight) }, .yellow);
        }
    }

    pub fn process(g: *Game) void {
        if (g.inputs.left) {
            g.player.vel[0] = -g.player.speedX;
            if (g.player.vel[0] < -g.player.maxvel[0]) g.player.vel[0] = -g.player.maxvel[0];
        } else if (g.inputs.right) {
            g.player.vel[0] = g.player.speedX;
            if (g.player.vel[0] > g.player.maxvel[0]) g.player.vel[0] = g.player.maxvel[0];
        } else {
            if (g.player.vel[0] > 0) g.player.vel[0] -= g.friction;
            if (g.player.vel[0] < 0) g.player.vel[0] += g.friction;
        }
        if (rl.isKeyPressed(rl.KeyboardKey.space)) {
            g.player.vel[1] = -g.player.jump_power;
        }

        const gr = if (g.player.vel[1] > 0) 2 * g.gravity * g.dt else g.gravity * g.dt;
        g.player.vel[1] = if (g.player.vel[1] <= g.player.maxvel[1]) g.player.vel[1] + gr else g.player.maxvel[1];
        T.movePlayer(g);
    }

    pub fn run(g: *Game) !void {
        rl.initWindow(g.screenWidth, g.screenHeight, "Platformer 2d");
        defer rl.closeWindow(); // Close window and OpenGL context
        rl.setExitKey(rl.KeyboardKey.null);
        rl.setTargetFPS(g.fps); // Set our game to run at 60 frames-per-second
        while (!rl.windowShouldClose()) { // Detect window close button or ESC key
            // g.collision_frame_id += 1;
            g.updateInputs();
            g._testInput();
            if (!g.pause) {
                // g.process();
                try g.vision.updatePlayerVision();
            }

            rl.beginDrawing();
            defer rl.endDrawing();
            rl.clearBackground(.black);
            g.draw();
        }
    }
    pub fn deinit(g: *Game) void {
        g.vision.deinit();
        g.wmap.deinit(g.allocator);
        g.allocator.destroy(g);
    }
};
