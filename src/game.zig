const rl = @import("raylib");
const std = @import("std");
const T = @import("types.zig");
const Vision = @import("vision.zig").Vision;
const math = std.math;
const print = std.debug.print;

pub const Game = struct {
    shader: rl.Shader = undefined,
    player_pos_loc: i32 = undefined,
    radius_loc: i32 = undefined,
    resolution_loc: i32 = undefined,
    renderTexture: rl.RenderTexture2D = undefined,

    collision_step_id: u64 = 0,
    vision: *Vision = undefined,
    pause: bool = false,
    fps: i32 = 60,
    gravity: f32 = 600,
    friction: f32 = 10,
    dt: f32 = undefined,
    allocator: std.mem.Allocator,
    player: T.Player = undefined,
    wmap: *T.WorldMap = undefined,
    screenWidth: i32 = T.fToI32(T.WorldMap.TileSize * T.WorldMap.TileNumberX),
    screenHeight: i32 = T.fToI32(T.WorldMap.TileSize * T.WorldMap.TileNumberY),

    inputs: struct {
        right: bool,
        left: bool,
        up: bool,
        down: bool,
        jump: bool,
        dash: bool,
        attack: bool,
        escape: bool,
        space: bool,
        space_release: bool,
    } = .{
        .right = false,
        .left = false,
        .up = false,
        .down = false,
        .jump = false,
        .dash = false,
        .attack = false,
        .escape = false,
        .space = false,
        .space_release = false,
    },

    pub fn init(allocator: std.mem.Allocator) !*Game {
        const game = try allocator.create(Game);
        game.* = .{
            .allocator = allocator,
        };
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
        g.inputs.space = rl.isKeyPressed(rl.KeyboardKey.space);
        g.inputs.space_release = rl.isKeyReleased(rl.KeyboardKey.space);
        g.inputs.escape = rl.isKeyPressed(rl.KeyboardKey.escape);
        if (g.inputs.escape) {
            g.pause = !g.pause;
            print("Pause: {}\n", .{g.pause});
        }
    }

    fn loadTextures(g:*Game) !void{
        g.player.texture_idle= try rl.loadTexture("assets/sprites/player-idle.png");
    }

    fn setup(g: *Game) !void {
        g.wmap = try T.WorldMap.init(g.allocator);
        g.vision = try Vision.init(g);
        g.dt = 1.0 / T.iToF32(g.fps);
        g.player = T.Player{
            .pos = .{ T.iToF32(@divTrunc(g.screenWidth, 2)), T.iToF32(@divTrunc(g.screenHeight, 2)) },
        };

        // Shader setup
        g.shader = try rl.loadShader(null, "assets/shaders/vision2.glsl");
        if (g.shader.id == 0) {
            print("ERROR: Failed to load shader\n", .{});
            // Handle error or return
        }

        // Get Uniform Locations
        g.player_pos_loc = rl.getShaderLocation(g.shader, "player_pos");
        g.radius_loc = rl.getShaderLocation(g.shader, "radius");
        g.resolution_loc = rl.getShaderLocation(g.shader, "resolution");

        // Set constant uniforms (Resolution doesn't change)
        const res = [2]f32{ T.iToF32(g.screenWidth), T.iToF32(g.screenHeight) };
        rl.setShaderValue(g.shader, g.resolution_loc, &res, rl.ShaderUniformDataType.vec2);

        // Load Render Texture (Off-screen canvas)
        g.renderTexture = try rl.loadRenderTexture(g.screenWidth, g.screenHeight);
        // rl.setTextureFilter(g.renderTexture.texture, rl.TextureFilter.bilinear);
        // try g.loadTextures();
    }

    fn draw(g: *Game) void {
        rl.beginTextureMode(g.renderTexture);
        rl.clearBackground(rl.Color.blank);
        g.vision.drawPlayerVision();
        rl.endTextureMode();
        if (g.pause) g.drawGridLines();
        for (g.wmap.platforms.items) |p| {
            if (!p.drawable) continue;
            const left = @max(p.pos[0] * T.WorldMap.TileSize, g.player.pos[0] * T.WorldMap.TileSize - g.player.vision_r * T.WorldMap.TileSize, 0);
            const right = @min(p.pos[0] * T.WorldMap.TileSize + p.size[0] * T.WorldMap.TileSize, g.player.pos[0] * T.WorldMap.TileSize + g.player.vision_r * T.WorldMap.TileSize, T.iToF32(g.screenWidth));
            const top = @max(p.pos[1] * T.WorldMap.TileSize, g.player.pos[1] * T.WorldMap.TileSize - g.player.vision_r * T.WorldMap.TileSize, 0);
            const bottom = @min(p.pos[1] * T.WorldMap.TileSize + p.size[1] * T.WorldMap.TileSize, g.player.pos[1] * T.WorldMap.TileSize + g.player.vision_r * T.WorldMap.TileSize, T.iToF32(g.screenHeight));
            if (left < right and top < bottom)
                rl.drawRectangleV(.{ .x = left, .y = top }, .{ .x = right - left, .y = bottom - top }, p.color);
        }
        g.player.draw();
        // g.player.drawTexture();
        rl.beginShaderMode(g.shader);
        const p_pos = [2]f32{ g.player.pos[0], g.player.pos[1] };
        const rad = g.player.vision_r;
        rl.setShaderValue(g.shader, g.player_pos_loc, &p_pos, rl.ShaderUniformDataType.vec2);
        rl.setShaderValue(g.shader, g.radius_loc, &rad, rl.ShaderUniformDataType.float);
        const tex = g.renderTexture.texture;
        rl.drawTextureRec(tex, rl.Rectangle{ .x = 0, .y = 0, .width = T.iToF32(tex.width), .height = -T.iToF32(tex.height) }, rl.Vector2{ .x = 0, .y = 0 }, rl.Color.white);
        rl.endShaderMode();
    }

    fn drawGridLines(g: *Game) void {
        var i: i32 = 0;
        while (i < T.WorldMap.TileNumberY) : (i += 1) {
            rl.drawLineV(rl.Vector2{ .x = 0, .y = T.iToF32(i) * T.WorldMap.TileSize }, rl.Vector2{ .x = T.iToF32(g.screenWidth), .y = T.iToF32(i) * T.WorldMap.TileSize }, .yellow);
        }
        i = 0;
        while (i < T.WorldMap.TileNumberX) : (i += 1) {
            rl.drawLineV(rl.Vector2{ .x = T.iToF32(i) * T.WorldMap.TileSize, .y = 0 }, rl.Vector2{ .x = T.iToF32(i) * T.WorldMap.TileSize, .y = T.iToF32(g.screenHeight) }, .yellow);
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
        if (g.inputs.space) {
            g.player.vel[1] = -g.player.jump_power;
        }
        if (g.inputs.space_release and g.player.vel[1] < 0) g.player.vel[1] /= 2;

        const gr = if (g.player.vel[1] > 0) 2 * g.gravity * g.dt else g.gravity * g.dt;
        g.player.vel[1] = if (g.player.vel[1] <= g.player.maxvel[1]) g.player.vel[1] + gr else g.player.maxvel[1];
        T.movePlayer(g);
    }

    pub fn run(g: *Game) !void {
        rl.initWindow(g.screenWidth, g.screenHeight, "Platformer 2d");
        defer rl.closeWindow(); // Close window and OpenGL context
        try g.setup();
        rl.setExitKey(rl.KeyboardKey.null);
        rl.setTargetFPS(g.fps); // Set our game to run at 60 frames-per-second
        while (!rl.windowShouldClose()) { // Detect window close button or ESC key
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
        rl.unloadRenderTexture(g.renderTexture);
        rl.unloadShader(g.shader);
        g.vision.deinit();
        g.wmap.deinit(g.allocator);
        g.allocator.destroy(g);
    }
};
