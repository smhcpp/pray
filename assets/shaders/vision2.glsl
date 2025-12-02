#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
out vec4 finalColor;
uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform vec2 player_pos;
uniform float radius;
uniform vec2 resolution;

void main()
{
    // 1. Get the Mask Value from the texture
    // We drew solid WHITE triangles in the texture.
    // Inside the triangle (visible area): .r is 1.0
    // Outside the triangle (hidden area): .r is 0.0
    float mask = texture(texture0, fragTexCoord).r;

    // 2. Calculate Distance for Soft Edges
    vec2 pixelPos = fragTexCoord * resolution;

    // Fix: Flip Y because OpenGL (0,0) is bottom-left, but Raylib/Player is top-left
    vec2 corrected_player_pos = vec2(player_pos.x, resolution.y - player_pos.y);

    float dist = distance(pixelPos, corrected_player_pos);

    float minimum = 0.4;
    // 3. Calculate Visibility Factor (0.0 to 1.0)
    // Falloff: 1.0 at center, fading to 0.0 at radius edge
    float falloff = clamp(1.0 - (dist / radius), minimum, 1.0);

    // Smooth the falloff curve (Hermite interpolation) for nicer looking light
    // falloff = falloff * falloff * (3.0 - 2.0 * falloff);
    if (mask < minimum) mask = minimum;

    // A pixel is visible ONLY if it is inside the mask AND within radius
    float visibility = mask * falloff;

    // 4. Output the SHADOW
    // Color is BLACK (0,0,0).
    // Alpha determines how dark the shadow is.
    // If visibility is 1.0 (High) -> Alpha is 0.0 (Transparent/No Shadow/See Game)
    // If visibility is 0.0 (Low)  -> Alpha is 1.0 (Opaque Black Shadow/Hidden)
    finalColor = vec4(0.0, 0.0, 0.0, 1.0 - visibility);
}
