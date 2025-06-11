package bcc.javaJostle;

import java.awt.Graphics;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Map {
    private int[][] tiles;

    public Map(String name) {
        // load the map from a file or resource
        ArrayList<int[]> rows = new ArrayList<>();
        String path = "app/src/main/resources/maps/" + name + ".txt"; // Ensure this path is correct
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.trim().split("\\s+");
                int[] row = new int[tokens.length];
                for (int i = 0; i < tokens.length; i++) {
                    switch (tokens[i]) {
                        case "W":
                            row[i] = Utilities.WALL;
                            break;
                        case "G":
                            row[i] = Utilities.GRASS;
                            break;
                        case "M":
                            row[i] = Utilities.MUD;
                            break;
                        default:
                            row[i] = -1; // unknown
                    }
                }
                rows.add(row);
            }
            tiles = rows.toArray(new int[rows.size()][]);
        } catch (IOException e) {
            e.printStackTrace();
            tiles = new int[0][0]; // Initialize to empty if loading fails
        }
    }

    public int[][] getTiles() {
        return tiles;
    }

    public void display(Graphics g, int panelWidth, int panelHeight, int cameraX, int cameraY, double zoomFactor, boolean[][] dangerMap) {
        if (tiles == null || tiles.length == 0) {
            return;
        }
        int rows = tiles.length;
        int cols = tiles[0].length;

        double currentTileSize = Utilities.TILE_SIZE * zoomFactor;

        // Calculate the range of tiles to draw based on camera and zoom
        int startCol = (int) Math.floor(cameraX / currentTileSize);
        int endCol = (int) Math.ceil((cameraX + panelWidth) / currentTileSize);
        int startRow = (int) Math.floor(cameraY / currentTileSize);
        int endRow = (int) Math.ceil((cameraY + panelHeight) / currentTileSize);

        // Clamp the ranges to the map boundaries
        startCol = Math.max(0, startCol);
        endCol = Math.min(cols, endCol);
        startRow = Math.max(0, startRow);
        endRow = Math.min(rows, endRow);

        for (int r = startRow; r < endRow; r++) {
            for (int c = startCol; c < endCol; c++) {
                int tileType = tiles[r][c];
                if (tileType != -1) {
                    // Calculate screen position of the tile
                    double x = (c * currentTileSize) - cameraX;
                    double y = (r * currentTileSize) - cameraY;

                    // No need for an additional check here as the loop bounds already handle culling.
                    // The previous check was:
                    // if (x + currentTileSize > 0 && x < panelWidth && y + currentTileSize > 0 && y < panelHeight) {
                    java.awt.Image img = null;
                    if(dangerMap != null && dangerMap.length > 1 && dangerMap[c][r]) {
                        tileType = Utilities.Danger;
                    }
                    switch (tileType) {
                        case Utilities.Danger: //danger (testing)
                            img = Utilities.DANGER_IMAGE;
                            break;
                        case Utilities.WALL:
                            img = Utilities.WALL_IMAGE;
                            break;
                        case Utilities.GRASS:
                            img = Utilities.GRASS_IMAGE;
                            break;
                        case Utilities.MUD:
                            img = Utilities.MUD_IMAGE;
                            break;
                    }
                    if (img != null) {
                        // Draw the image at the calculated screen position
                        g.drawImage(img, (int) x, (int) y, (int) currentTileSize, (int) currentTileSize, null);
                    }
                    // }
                }
            }
        }
    }
}