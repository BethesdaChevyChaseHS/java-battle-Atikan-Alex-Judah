package bcc.javaJostle;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;

import javax.imageio.ImageIO;

public class Utilities {
    public static BufferedImage DANGER_IMAGE;
    public static BufferedImage WALL_IMAGE;
    public static BufferedImage GRASS_IMAGE;
    public static BufferedImage MUD_IMAGE;
    public static BufferedImage ROBOT_ERROR;
    public static BufferedImage DEFAULT_PROJECTILE_IMAGE;
    public static BufferedImage HEALTH_PACK_IMAGE;
    public static BufferedImage SPEED_PACK_IMAGE;
    public static BufferedImage ATTACK_PACK_IMAGE;

    public final static int WALL = 0;
    public final static int GRASS = 1;
    public final static int MUD = 2;
     public final static int Danger = 3;

    public final static double POWER_UP_SPAWN_CHANCE = .003;

    public static final int SCREEN_WIDTH = 800;
    public static final int SCREEN_HEIGHT = 600;
    public static final int TILE_SIZE = 32; // Default tile size
    public static final int PROJECTILE_SIZE = 10; // Default projectile size
    public static final int POWER_UP_SIZE = 20; // Default power-up size
    public static final int ROBOT_SIZE = 28; // Default robot size

    public static ArrayList<Integer> keysPressed = new ArrayList<Integer>();

    public static BufferedImage loadImage(String imgName) {
        try {
            // Ensure the full path is constructed correctly
            BufferedImage img = ImageIO.read(new File("app/src/main/resources/images/" + imgName));
            return cropToContent(img);
        } catch (IOException e) {
            System.err.println("Failed to load image: " + imgName + " - " + e.getMessage());
            e.printStackTrace();
            return null; // Return null if loading fails
        }
    }

    public static void loadImages() {
        // Use the loadImage method for each image
        DANGER_IMAGE = loadImage("dangerMap.png");
        WALL_IMAGE = loadImage("wall.png");
        GRASS_IMAGE = loadImage("grass.png");
        MUD_IMAGE = loadImage("mud.png");
        ROBOT_ERROR = loadImage("robotError.png");
        DEFAULT_PROJECTILE_IMAGE = loadImage("defaultProjectile.png");

        // Optional: Add checks here if any image failed to load, though loadImage
        // already prints errors
        if (WALL_IMAGE == null) {
            System.err.println("WALL_IMAGE could not be loaded.");
        }
        if (GRASS_IMAGE == null) {
            System.err.println("GRASS_IMAGE could not be loaded.");
        }
        if (MUD_IMAGE == null) {
            System.err.println("MUD_IMAGE could not be loaded.");
        }
    }

    public static BufferedImage cropToContent(BufferedImage src) {
        int width = src.getWidth();
        int height = src.getHeight();

        int minX = width;
        int minY = height;
        int maxX = 0;
        int maxY = 0;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int pixel = src.getRGB(x, y);
                int alpha = (pixel >> 24) & 0xff;

                if (alpha > 0) {
                    if (x < minX)
                        minX = x;
                    if (y < minY)
                        minY = y;
                    if (x > maxX)
                        maxX = x;
                    if (y > maxY)
                        maxY = y;
                }
            }
        }

        if (maxX < minX || maxY < minY) {
            // No visible content
            return new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB);
        }

        return src.getSubimage(minX, minY, (maxX - minX + 1), (maxY - minY + 1));
    }

    public static void handleKeyPressed(int keyCode) {
        if (!keysPressed.contains(keyCode)) {
            keysPressed.add(keyCode);
        }
    }

    public static void handleKeyReleased(int keyCode) {
        for (int i = 0; i < keysPressed.size(); i++) {
            if (keysPressed.get(i) == keyCode) {
                keysPressed.remove(i); // Remove using index for clarity
                break;
            }
        }
    }

    // For testing purposes only
    private static boolean[] testKeyStates = new boolean[256];

    // Test helper methods
    public static void resetKeyStates() {
        for (int i = 0; i < testKeyStates.length; i++) {
            testKeyStates[i] = false;
        }
    }

    public static void setKeyPressed(int keyCode, boolean isPressed) {
        if (keyCode < testKeyStates.length) {
            testKeyStates[keyCode] = isPressed;
        }
    }

    // Modify the existing isKeyPressed method to use the test states in test mode
    private static boolean inTestMode = false;

    public static void setTestMode(boolean enabled) {
        inTestMode = enabled;
    }

    // Update the existing isKeyPressed method
    public static boolean isKeyPressed(int keyCode) {
        if (inTestMode) {
            return keyCode < testKeyStates.length && testKeyStates[keyCode];
        } else {
            // Original implementation
            return keysPressed.contains(keyCode);
        }
    }

    public static Robot createRobot(int x, int y, String className) {
        // The URL for the '/robots' resource directory.
        // This will be used if className is not "MyRobot".
       
        // Attempt to define robotsResourceRootUrl from a specific source directory
        File specificSrcRobotsDir = new File("app/src/main/resources/robots"); // Path relative to project root
        URL robotsResourceRootUrl = null;


        if (specificSrcRobotsDir.exists() && specificSrcRobotsDir.isDirectory()) {
            try {
                robotsResourceRootUrl = specificSrcRobotsDir.toURI().toURL();
                // Optional: Add a log to confirm which path is being used
                // System.out.println("Attempting to use direct src path for robots: " + robotsResourceRootUrl);
            } catch (java.net.MalformedURLException e) {
                System.err.println("Error creating URL for src robots directory: " + specificSrcRobotsDir.getAbsolutePath() + " - " + e.getMessage());
                // robotsResourceRootUrl will remain null, handled by later checks
            }
        }

        if (className.equals("MyRobot")) {
            // Assuming MyRobot is always on the standard classpath
            return new MyRobot(x, y);
        } else {
            // For other robot classes, attempt to load them using URLClassLoader
            // from the '/robots' resource directory.
            if (robotsResourceRootUrl == null) {
                System.err.println("Robots resource directory '/robots' not found in classpath. Cannot dynamically load: " + className);
                // As a fallback, try to load the class directly from the classpath.
                // This might be useful if some robots are packaged normally.
                try {
                    System.out.println("Fallback: Attempting to load robot class directly from classpath: bcc.javaJostle." + className);
                    Class<?> loadedClass = Class.forName("bcc.javaJostle." + className);
                    if (Robot.class.isAssignableFrom(loadedClass)) {
                        Constructor<?> constructor = loadedClass.getConstructor(int.class, int.class);
                        return (Robot) constructor.newInstance(x, y);
                    } else {
                        System.err.println("Class bcc.javaJostle." + className + " (loaded from classpath) does not extend Robot.");
                        return null; // Class found but not a Robot or doesn't have the right constructor
                    }
                } catch (Exception directLoadException) {
                    System.err.println("Fallback failed: Could not load bcc.javaJostle." + className + " directly from classpath: " + directLoadException.getMessage());
                    return null; // Fallback loading failed
                }
            }


            // If robotsResourceRootUrl is valid, proceed with URLClassLoader
            try {
                URL[] urls = { robotsResourceRootUrl };
                // Using Utilities.class.getClassLoader() as the parent.
                // The original code used Game.class.getClassLoader().
                URLClassLoader classLoader = new URLClassLoader(urls, Utilities.class.getClassLoader());
               
                System.out.println("Attempting to load robot class via URLClassLoader: bcc.javaJostle." + className + " from " + robotsResourceRootUrl);
                // Assumes robot classes are in the 'bcc.javaJostle' package within the /robots directory
                Class<?> loadedClass = classLoader.loadClass("bcc.javaJostle." + className);


                if (Robot.class.isAssignableFrom(loadedClass)) {
                    Constructor<?> constructor = loadedClass.getConstructor(int.class, int.class);
                    Robot robot = (Robot) constructor.newInstance(x, y);
                    return robot;
                } else {
                    System.err.println("Class " + className + " (loaded via URLClassLoader from " + robotsResourceRootUrl + ") does not extend Robot.");
                }
            } catch (ClassNotFoundException e) {
                System.err.println("Robot class not found via URLClassLoader: bcc.javaJostle." + className + " in " + robotsResourceRootUrl + " - " + e.getMessage());
            } catch (NoSuchMethodException e) {
                System.err.println("Constructor (int, int) not found for " + className + " (loaded via URLClassLoader from " + robotsResourceRootUrl + ") - " + e.getMessage());
            } catch (Exception e) {
                System.err.println("Error loading or instantiating " + className + " via URLClassLoader from " + robotsResourceRootUrl + ": " + e.getMessage());
                e.printStackTrace();
            }
        }
        // Return null if MyRobot was not the class, and dynamic/fallback loading failed.
        return null;
    }

}