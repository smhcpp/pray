const std = @import("std");
const print = std.debug.print;
const Game = @import("game.zig").Game;
pub fn main() anyerror!void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();
    const g = try Game.init(allocator);
    try g.run();
    g.deinit();
    const leak = gpa.deinit(); // Checks for leaks in debug
    print("\nLeaks:\n {}", .{leak});
}
