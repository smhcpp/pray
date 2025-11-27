#version 330

in vec2 fragTexCoord;
in vec4 fragColor;

out vec4 finalColor;

uniform sampler2D texture0;
uniform vec2 player_pos;
uniform float radius;
uniform vec2 resolution;

void main() {
    vec4 texColor = texture(texture0, fragTexCoord);

    // Convert to pixel coordinates
    vec2 pixelPos = fragTexCoord * resolution;

    float dist = distance(pixelPos, player_pos);

    if (dist > radius) {
        finalColor = vec4(0.0, 0.0, 0.0, 1.0); // Black outside
    } else {
        finalColor = texColor; // Normal inside
    }
}
