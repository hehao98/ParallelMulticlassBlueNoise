

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>

// OpenGL Math Library
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "stb_image.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw_gl3.h"

#include "Shader.h"
#include "PointSet.h"
#include "SpectrumImage.h"

// **********GLFW window related functions**********
// Returns pointer to a initialized window with OpenGL context set up
GLFWwindow *init();
// Sometimes user might resize the window. so the OpenGL viewport should be adjusted as well.
void frameBufferSizeCallback(GLFWwindow *window, int width, int height);
// User input is handled in this function
void processInput(GLFWwindow *window);
// Mouse input is handled in this function
void mouseCallback(GLFWwindow *window, double xpos, double ypos);
void scrollCallback(GLFWwindow *window, double offsetX, double offsetY);
// ********** ImGui Utilities **********
void imGuiInit(GLFWwindow *window);
void imGuiSetup(GLFWwindow *window);

int gScreenWidth = 1200;
int gScreenHeight = 600;
float gDeltaTime = 0.0f;
float gLastFrame = 0.0f;
int classType = -1;

int main()
{
    GLFWwindow *window = init();
    if (window == nullptr) {
        std::cout << "Failed to initialize GLFW and OpenGL!" << std::endl;
        return -1;
    }

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    imGuiInit(window);

    PointSet pointSet;
    //pointSet.generateWhiteNoisePointSet(10000, 0, 500, 0, 500);
    pointSet.readPointSetFromFile("../testdata/pointset1.txt", 0, 1, 0, 1);
    pointSet.updateRenderData(-1, 0, -1, 1);

    SpectrumImage image;
    image.readSpectrumFromFile("../testdata/spectrum2.txt");
    //image.updateSpectrum(pointSet.points);

    Shader shader("../shaders/point.vert", "../shaders/point.frag");
    Shader imageShader("../shaders/image.vert", "../shaders/image.frag");

    // Game loop
    while (!glfwWindowShouldClose(window)) {
        // Calculate how much time since last frame
        auto currentFrame = (float)glfwGetTime();
        gDeltaTime = currentFrame - gLastFrame;
        gLastFrame = currentFrame;

        processInput(window);

        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        imGuiSetup(window);

        pointSet.updateRenderData(-1, 0, -1, 1, classType);
        pointSet.render(shader);

        image.render(imageShader);

        // Render GUI last
        ImGui::Render();
        ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();

        // An ugly workaround to make rendering work on Mac OS Mojave
        int xpos, ypos;
        glfwGetWindowPos(window, &xpos, &ypos);
        glfwSetWindowPos(window, xpos+1, ypos+1);
        glfwSetWindowPos(window, xpos, ypos);
    }

    glfwTerminate();
    return 0;
}

void imGuiInit(GLFWwindow *window)
{
    // Setup ImGui binding
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable Gamepad Controls
	ImGui_ImplGlfwGL3_Init(window, true);
	// Setup style
	ImGui::StyleColorsDark();
}

void imGuiSetup(GLFWwindow *window)
{
    ImGui_ImplGlfwGL3_NewFrame();

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)",
                1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

    ImGui::InputInt("Enter Class Type", &classType);

    if (ImGui::Button("Quit")) {
        glfwSetWindowShouldClose(window, true);
    }
}

GLFWwindow *init()
{
    // Initialization of GLFW context
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // Something needed for Mac OS X
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    // Create a window object
    GLFWwindow *window = glfwCreateWindow(gScreenWidth, gScreenHeight, "Window Title", nullptr, nullptr);
    if (window == nullptr) {
        std::cout << "Failed to create GLFW window!" << std::endl;
        glfwTerminate();
        return nullptr;
    }

    glfwMakeContextCurrent(window);

    // Initialize GLAD before calling OpenGL functions
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return nullptr;
    }

    // Tell OpenGL the size of rendering window
    glViewport(0, 0, gScreenWidth, gScreenHeight);

    // Set the windows resize callback function
    glfwSetFramebufferSizeCallback(window, frameBufferSizeCallback);

    // Set up mouse input
    glfwSetCursorPosCallback(window, mouseCallback);
    glfwSetScrollCallback(window, scrollCallback);

    return window;
}

void frameBufferSizeCallback(GLFWwindow *window, int width, int height)
{
    gScreenWidth = width;
    gScreenHeight = height;
    glViewport(0, 0, gScreenWidth, gScreenHeight);
}

void processInput(GLFWwindow *window)
{
    // add this to avoid repetitive pressing
    static double timeLastPressed = 0.0;
    if (glfwGetTime() - timeLastPressed < 0.1)
        return;
    timeLastPressed = glfwGetTime();

    // Exit
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        timeLastPressed = glfwGetTime();
    }
}

void mouseCallback(GLFWwindow *window, double xpos, double ypos)
{

}

void scrollCallback(GLFWwindow *window, double offsetX, double offsetY)
{

}
