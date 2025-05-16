#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <GLFW/glfw3native.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#define NOMINMAX
#include <windows.h>

const int N = 8;
const float cellSize = 100.0f;
const float epsilon = 1e-3f;

glm::vec2 velocity[N][N];
bool showVelocity = true;

std::vector<float> density(N *N);
std::vector<float> prevDensity(N *N);
std::vector<float> densityVertices;

int windowWidth = static_cast<int>(N * cellSize);
int windowHeight = static_cast<int>(N * cellSize);

std::vector<glm::vec2> lineVertices;
GLuint lineVAO, lineVBO;
GLuint densityVAO, densityVBO;
GLuint lineShaderProgram, densityShaderProgram;
glm::mat4 projection = glm::ortho(0.0f, (float)windowWidth, 0.0f, (float)windowHeight);

inline int IX(int x, int y) { return x + y * N; }

enum BoundaryType
{
    BND_NONE = 0,
    BND_VELOCITY_X = 1,
    BND_VELOCITY_Y = 2,
    BND_DENSITY = 3
};

void setBoundary(int b, std::vector<float> &x)
{
    for (int i = 1; i < N - 1; ++i)
    {
        x[IX(0, i)] = (b == BND_VELOCITY_X) ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N - 1, i)] = (b == BND_VELOCITY_X) ? -x[IX(N - 2, i)] : x[IX(N - 2, i)];

        x[IX(i, 0)] = (b == BND_VELOCITY_Y) ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = (b == BND_VELOCITY_Y) ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

void initVelocityField()
{
    srand(static_cast<unsigned>(time(nullptr) + rand()));
    float angle = static_cast<float>(rand()) / RAND_MAX * 2.0f * 3.14159f;
    float speed = 0.1f;
    glm::vec2 dir = glm::vec2(cos(angle), sin(angle)) * speed;
    for (int x = 0; x < N; ++x)
    {
        for (int y = 0; y < N; ++y)
        {
            velocity[x][y] = dir;
        }
    }
}

void buildLineVertices()
{
    const float sizeMultiplier = cellSize * 2;
    lineVertices.clear();
    for (int x = 0; x < N; ++x)
    {
        for (int y = 0; y < N; ++y)
        {
            glm::vec2 base = glm::vec2(x, y) * cellSize + glm::vec2(cellSize / 2, cellSize / 2);
            glm::vec2 tip = base + velocity[x][y] * sizeMultiplier;

            lineVertices.push_back(base);
            lineVertices.push_back(tip);

            glm::vec2 dir = glm::normalize(velocity[x][y]);
            glm::vec2 ortho = glm::vec2(-dir.y, dir.x);

            float arrowSize = 3.0f;
            glm::vec2 left = tip - dir * arrowSize + ortho * arrowSize * 0.5f;
            glm::vec2 right = tip - dir * arrowSize - ortho * arrowSize * 0.5f;

            lineVertices.push_back(tip);
            lineVertices.push_back(left);

            lineVertices.push_back(tip);
            lineVertices.push_back(right);

            if (x < N - 1)
            { // Vertical grid lines (not at the last column)
                glm::vec2 start = glm::vec2((x + 1) * cellSize, 0.0f);
                glm::vec2 end = glm::vec2((x + 1) * cellSize, N * cellSize);
                lineVertices.push_back(start);
                lineVertices.push_back(end);
            }

            if (y < N - 1)
            { // Horizontal grid lines (not at the last row)
                glm::vec2 start = glm::vec2(0.0f, (y + 1) * cellSize);
                glm::vec2 end = glm::vec2(N * cellSize, (y + 1) * cellSize);
                lineVertices.push_back(start);
                lineVertices.push_back(end);
            }
        }
    }
}

void buildDensityVertices(const std::vector<float> &d)
{
    densityVertices.clear();
    for (int x = 0; x < N; ++x)
    {
        for (int y = 0; y < N; ++y)
        {
            float val = glm::clamp(d[IX(x, y)] / 100.0f, 0.0f, 1.0f);

            if (val < epsilon)
                continue;

            float px = x * cellSize + cellSize * 0.5f;
            float py = y * cellSize + cellSize * 0.5f;

            densityVertices.push_back(px);
            densityVertices.push_back(py);
            densityVertices.push_back(val);
            // std::cout << "(" << x << "," << y << "): "<< val << std::endl;
        }
    }

    glBindVertexArray(densityVAO);
    glBindBuffer(GL_ARRAY_BUFFER, densityVBO);
    glBufferData(GL_ARRAY_BUFFER, densityVertices.size() * sizeof(float), densityVertices.data(), GL_DYNAMIC_DRAW);
    glBindVertexArray(0);
}

void add_source(float *x, float *s, float dt)
{
    for (int i = 1; i < N * N; ++i)
    {
        x[i] += dt * s[i];
    }
}

// void diffuse(float* x, float* x0, float diff, float dt) {
//     float a = dt * diff * N * N;
//     int RELAXATIONS = 20;
//     for (int k = 0; k < RELAXATIONS; ++k)
//     {
//         for (int i = 1; i < N - 1; ++i)
//         {
//             for (int j = 1; j < N - 1; ++j)
//             {
//                 float divisor = 2 + 4 * a;
//                 x[IX(i,j)] = (x0[IX(i,j)] + a * gridSum(x, i, j)) / divisor;
//             }
//         }

//     }
// }

// void project(glm::vec2* vel, glm::vec2* velOut) {
//     std::vector<float> div(N * N, 0.0f);
//     std::vector<float> p(N * N, 0.0f);

//     float h = 1.0f / N;

//     for (int i = 1; i < N - 1; ++i)
//     {
//         for (int j = 1; j < N - 1; ++j)
//         {
//             float dx = vel[IX(i+1, j)].x - vel[IX(i-1,j)].x;
//             float dy = vel[IX(i,j+1)].y - vel[IX(i,j-1)].y;
//             div[IX(i,j)] = -0.5f * h * (dx + dy);
//             p[IX(i,j)] = 0.0f;
//         }
//     }

//     for (int i = 1; i < N - 1; ++i)
//     {
//         for (int j = 1; j < N -1; ++j)
//         {
//             p[IX(i,j)] = (div[IX(i,j)] + gridSum(p,i,j)) / 4.0f;
//         }

//     }

//     for (int i = 1; i < N-1; ++i)
//     {
//         for (int j = 1; j < N - 1; ++j)
//         {
//             float dx = vel[IX(i+1, j)].x - vel[IX(i-1,j)].x;
//             float dy = vel[IX(i,j+1)].y - vel[IX(i,j-1)].y;
//             velOut[IX(i,j)] = vel[IX(i,j)] - glm::vec2(dx,dy);
//         }

//     }
// }

float clampf(float x, float minVal, float maxVal)
{
    return std::max(minVal, std::min(x, maxVal));
}

void advect(std::vector<float> &d, std::vector<float> &d0, glm::vec2 velocity[N][N], float dt)
{
    float dt0 = dt * N;
    for (int i = 1; i < N - 1; ++i)
    {
        for (int j = 1; j < N - 1; ++j)
        {

            // Source location for the point at (i,j)
            glm::vec2 pos = glm::vec2(i, j) - dt0 * velocity[i][j];

            pos.x = clampf(pos.x, 0.5f, N - 1.5f);
            pos.y = clampf(pos.y, 0.5f, N - 1.5f);

            int i0 = (int)pos.x;
            int j0 = (int)pos.y;
            int i1 = i0 + 1;
            int j1 = j0 + 1;

            float sx = pos.x - i0;
            float sy = pos.y - j0;

            float d00 = d0[IX(i0, j0)];
            float d10 = d0[IX(i1, j0)];
            float d01 = d0[IX(i0, j1)];
            float d11 = d0[IX(i1, j1)];

            d[IX(i, j)] = (1 - sx) * ((1 - sy) * d00 + sy * d01) +
                          sx * ((1 - sy) * d10 + sy * d11);
        }
    }
}

void seedDensity(std::vector<float> &d)
{
    std::fill(d.begin(), d.end(), 0.0f);
    for (int i = 0; i < d.size() / 4; ++i)
    {
        int x = rand() % (N / 2) + N / 4;
        int y = rand() % (N / 2) + N / 4;
        d[IX(x, y)] = 100.0f;
    }
}

void clearDensity(std::vector<float> &d)
{
    std::fill(d.begin(), d.end(), 0.0f);
}

float gridSum(std::vector<float> vec, int i, int j)
{
    return vec[IX(i - 1, j)] + vec[IX(i + 1, j)] + vec[IX(i, j - 1)] + vec[IX(i, j + 1)];
}

void setupBuffers()
{
    // Line Vertices
    glGenVertexArrays(1, &lineVAO);
    glGenBuffers(1, &lineVBO);

    glBindVertexArray(lineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);

    glBufferData(GL_ARRAY_BUFFER, lineVertices.size() * sizeof(glm::vec2), lineVertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void *)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    // Density Vertices
    glGenVertexArrays(1, &densityVAO);
    glGenBuffers(1, &densityVBO);

    glBindVertexArray(densityVAO);
    glBindBuffer(GL_ARRAY_BUFFER, densityVBO);

    glBufferData(GL_ARRAY_BUFFER, densityVertices.size() * sizeof(float), densityVertices.data(), GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}

void renderFrame()
{
    glClearColor(0.1f, 0.1f, 0.1f, 0.1f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUseProgram(densityShaderProgram);
    glBindVertexArray(densityVAO);
    glUniformMatrix4fv(glGetUniformLocation(densityShaderProgram, "uProjection"), 1, GL_FALSE, glm::value_ptr(projection));

    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(densityVertices.size() / 3));
    glBindVertexArray(0);

    if (showVelocity)
    {
        glUseProgram(lineShaderProgram);
        glBindVertexArray(lineVAO);
        glUniformMatrix4fv(glGetUniformLocation(lineShaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(lineVertices.size()));
        glBindVertexArray(0);
    }
}

std::string loadFile(const std::string &path)
{
    std::ifstream file(path);
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

GLuint compileShader(GLenum type, const std::string &source)
{
    GLuint shader = glCreateShader(type);
    const char *src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        std::cerr << "Shader compilation error: " << log << std::endl;
    }

    return shader;
}

GLuint createShaderProgram(const std::string &vertPath, const std::string &fragPath)
{
    GLuint vert = compileShader(GL_VERTEX_SHADER, loadFile(vertPath));
    GLuint frag = compileShader(GL_FRAGMENT_SHADER, loadFile(fragPath));

    GLuint program = glCreateProgram();
    glAttachShader(program, vert);
    glAttachShader(program, frag);
    glLinkProgram(program);

    glDeleteShader(vert);
    glDeleteShader(frag);

    return program;
}

void keyCallBack(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS)
    {
        if (key == GLFW_KEY_V)
        {
            showVelocity = !showVelocity;
        }
    }
}

void mouseButtonCallBack(GLFWwindow *window, int button, int action, int mods)
{
    std::cout << "In mouse call back" << std::endl;
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        int i = static_cast<int>(xpos / cellSize);
        int j = static_cast<int>((windowHeight - ypos) / cellSize); // invert since origin is bottom-left
        std::cout << "If condition met cell: " << i << ", " << j << std::endl;
        if (i >= 0 && i < N && j >= 0 && j < N)
        {
            density[IX(i, j)] += 10.0f;
            prevDensity[IX(i, j)] += 10.0f; // Keep buffers in sync
        }
    }
}

void updateSimulation(float dt)
{
    std::swap(density, prevDensity);
    advect(density, prevDensity, velocity, dt);
    setBoundary(BND_DENSITY, density);
}

int main()
{

    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow *window = glfwCreateWindow(windowWidth, windowHeight, "Fluid Velocity Field", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    glfwSetKeyCallback(window, keyCallBack);
    glfwSetMouseButtonCallback(window, mouseButtonCallBack);
    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }
    glViewport(0, 0, windowWidth, windowHeight);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    lineShaderProgram = createShaderProgram("../../src/shaders/vertex.glsl", "../../src/shaders/fragment.glsl");
    densityShaderProgram = createShaderProgram("../../src/shaders/densityVertex.glsl", "../../src/shaders/densityFragment.glsl");

    initVelocityField();
    buildLineVertices();
    seedDensity(density);
    clearDensity(prevDensity);
    setupBuffers();

    float currTime = static_cast<float>(glfwGetTime());
    float lastTime = currTime;
    while (!glfwWindowShouldClose(window))
    {

        float newTime = static_cast<float>(glfwGetTime());
        float dt = newTime - lastTime;
        lastTime = newTime;
        // input processing

        // rendering
        updateSimulation(dt);
        buildDensityVertices(density);
        renderFrame();

        // events and buffer swap
        glfwSwapBuffers(window);
        glfwPollEvents();

        GLenum err;
        while ((err = glGetError()) != GL_NO_ERROR)
        {
            std::cout << "GL Error: " << err << std::endl;
        }

        // Sleep(1000);
    }
    glfwTerminate();
    return 0;
}