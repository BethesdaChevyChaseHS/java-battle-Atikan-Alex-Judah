package bcc.javaJostle;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

public class MyRobot extends Robot{
    private int lastX;
    private int lastY;
    private int distancesDist;
    private int lastBotX = -1;
    private int lastBotY = -1;
    private int targetTileX;
    private int targetTileY;
    private int[][] tiles; 
    private int width;
    private int height;
    private boolean[][] visited;
    private int[][] distances;
    public boolean[][] dangerMap;
    private boolean flip;
    private int stuckCounter = 0;

    public MyRobot(int x, int y){
        super(x, y, 2, 3, 2, 3,"bob", "myRobotImage.png", "defaultProjectileImage.png");
        // Health: 3, Speed: 3, Attack Speed: 2, Projectile Strength: 2
        // Total = 10
        // Image name is "myRobotImage.png"
    }

    private class Tile {
        private int x;
        private int y;

        public Tile(int x, int y) {
            this.x = x;
            this.y = y;
        }

        public int tileX() {
            return x;
        }

        public int tileY() {
            return y;
        }
    }

    public void think(ArrayList<Robot> robots, ArrayList<Projectile> projectiles, Map map, ArrayList<PowerUp> powerups) {
        int high = 0;
        ArrayList<Robot> otherRobots = new ArrayList<Robot>();
        for(int i = 0; i < robots.size(); i++) {
            if(robots.get(i) != this && robots.get(i).isAlive()) {
                if(robots.get(i).getProjectileStrengthPoints() + 5 > high) { 
                    high = robots.get(i).getProjectileStrengthPoints() + 5;
                }
                otherRobots.add(robots.get(i));
            }
        }
        distancesDist = high - getSpeed();
        if(powerups.size() > 0) {
            dijkstras(map, projectiles, powerups.get(0).getX(), powerups.get(0).getY(), otherRobots);
        } else {
            farDijkstras(otherRobots, projectiles, map);
        }
        if(canAttack()){
            for(Robot robot : robots) {
                if (robot != this && robot.isAlive() ){
                    if(predictiveAimAndShoot(robot)) {
                        break;
                    }
                }
            }
        }
        int currentX = getX();
        int currentY = getY();
        if (currentX == lastBotX && currentY == lastBotY) {
            stuckCounter++;
        } else {
            stuckCounter = 0;
        }
        lastBotX = currentX;
        lastBotY = currentY;
        if (stuckCounter >= 2) {
            int[][] tiles = map.getTiles();
            int width = tiles[0].length;
            int height = tiles.length;
            for (int tries = 0; tries < 10; tries++) {
                int rx = (int)(Math.random() * width);
                int ry = (int)(Math.random() * height);
                if (tiles[ry][rx] != Utilities.WALL) {
                    dijkstras(map, projectiles, rx * Utilities.TILE_SIZE, ry * Utilities.TILE_SIZE, otherRobots);
                    break;
                }
            }
            stuckCounter = 0;
        }
    }

    private void farDijkstras(ArrayList<Robot> robots, ArrayList<Projectile> projectiles, Map map) {
        ArrayList<Robot> otherRobots = new ArrayList<Robot>();
        tiles = map.getTiles();
        width = map.getTiles()[0].length;
        height = map.getTiles().length;
        visited = new boolean[width][height];
        distances = new int[width][height];

        for(Robot robot : robots) {
            if(robot != this) {
                otherRobots.add(robot);
            }
        }

        for(Robot robot : otherRobots) {
            for(int i = 0; i < width; i++) {
                for(int j = 0; j < height; j++) {
                    if(tiles[j][i] == Utilities.WALL) {  
                        visited[i][j] = true;
                    }else{
                        visited[i][j] = false;
                    }
                    if(robot.getX()/Utilities.TILE_SIZE == i && robot.getY()/Utilities.TILE_SIZE == j) {
                        distances[i][j] = 0;
                    }else{
                        distances[i][j] = Integer.MAX_VALUE;
                    }
                }
            }

            Queue<Tile> priority = new LinkedList<>();
            priority.add(new Tile(robot.getX()/Utilities.TILE_SIZE, robot.getY()/Utilities.TILE_SIZE));
            visited[robot.getX()/Utilities.TILE_SIZE][robot.getY()/Utilities.TILE_SIZE] = true;
            while(!priority.isEmpty()) {
                int curr = priority.size();
                for(int i = 0; i < curr; i++) {
                    Tile tile = priority.poll();
                    int tileX = tile.tileX();
                    int tileY = tile.tileY();
                    if(tileX + 1 < visited.length && !visited[tileX + 1][tileY]) {
                        visited[tileX + 1][tileY] = true;
                        distances[tileX + 1][tileY] = Math.min(distances[tileX][tileY] + (tiles[tileY][tileX + 1] == Utilities.MUD ? 2 : 1), distances[tileX + 1][tileY]);
                        priority.add(new Tile(tileX + 1, tileY));
                    }
                    if(tileX - 1 >= 0 && !visited[tileX - 1][tileY]) {
                        visited[tileX - 1][tileY] = true;
                        distances[tileX - 1][tileY] = Math.min(distances[tileX][tileY] + (tiles[tileY][tileX - 1] == Utilities.MUD ? 2 : 1), distances[tileX - 1][tileY]);
                        priority.add(new Tile(tileX - 1, tileY));
                    }
                    if(tileY + 1 < visited[0].length && !visited[tileX][tileY + 1]) {
                        visited[tileX][tileY + 1] = true;
                        distances[tileX][tileY + 1] = Math.min(distances[tileX][tileY] + (tiles[tileY + 1][tileX] == Utilities.MUD ? 2 : 1), distances[tileX][tileY + 1]);
                        priority.add(new Tile(tileX, tileY + 1));
                    }
                    if(tileY - 1 >= 0 && !visited[tileX][tileY - 1]) {
                        visited[tileX][tileY - 1] = true;
                        distances[tileX][tileY - 1] = Math.min(distances[tileX][tileY] + (tiles[tileY - 1][tileX] == Utilities.MUD ? 2 : 1), distances[tileX][tileY - 1]);
                        priority.add(new Tile(tileX, tileY - 1));
                    }
                }
            }
        }

        int farthestX = 0;
        int farthestY = 0;
        int farthestDist = Integer.MIN_VALUE;

        for(int i = 1; i < distances.length - 1; i++) {
            for(int j = 1; j < distances[0].length - 1; j++) {
                if(distances[i][j] > farthestDist && distances[i][j] != Integer.MAX_VALUE) {
                    farthestDist = distances[i][j];
                    farthestX = i;
                    farthestY = j;
                }
            }
        }
        
        dijkstras(map, projectiles, (farthestX + (flip ? (farthestX < 12 ? 2 : -2) : 0)) * Utilities.TILE_SIZE, (farthestY + (flip ? (farthestY < 9 ? 2 : -2) : 0)) * Utilities.TILE_SIZE, otherRobots);
    } 

    private boolean predictiveAimAndShoot(Robot target) {
        if(!hasLineOfSight(target)) {
            return false;
        }
        double shooterX = getX();
        double shooterY = getY();
        double targetX = target.getX();
        double targetY = target.getY();
        double targetVelX = target.xMovement * target.getSpeed();
        double targetVelY = target.yMovement * target.getSpeed();
        double projectileSpeed = getProjectileStrengthPoints() + 5;
        double relX = targetX - shooterX;
        double relY = targetY - shooterY;
        double a = targetVelX * targetVelX + targetVelY * targetVelY - projectileSpeed * projectileSpeed;
        double b = 2 * (relX * targetVelX + relY * targetVelY);
        double c = relX * relX + relY * relY;
        double discriminant = b * b - 4 * a * c;
        double interceptTime;
        if (a == 0) {
            interceptTime = -c / b;
        } else if (discriminant >= 0) {
            double sqrtDisc = Math.sqrt(discriminant);
            double t1 = (-b + sqrtDisc) / (2 * a);
            double t2 = (-b - sqrtDisc) / (2 * a);
            interceptTime = Math.min(t1, t2);
            if (interceptTime < 0) interceptTime = Math.max(t1, t2);
        } else {
            interceptTime = 0;
        }
        if (interceptTime < 0) interceptTime = 0;
        double aimX = targetX + targetVelX * interceptTime;
        double aimY = targetY + targetVelY * interceptTime;
        shootAtLocation((int) aimX, (int) aimY);
        return true;
    }

    public void dijkstras(Map map, ArrayList<Projectile> projectiles, double targetX, double targetY, ArrayList<Robot> otherRobots) {
        tiles = map.getTiles();
        width = map.getTiles()[0].length;
        height = map.getTiles().length;
        visited = new boolean[width][height];
        distances = new int[width][height];

        int robotTileX = getX() / Utilities.TILE_SIZE;
        int robotTileY = getY() / Utilities.TILE_SIZE;
        int robotTileXPlus = (getX() + Utilities.ROBOT_SIZE) / Utilities.TILE_SIZE;
        int robotTileYPlus = (getY() + Utilities.ROBOT_SIZE) / Utilities.TILE_SIZE;
        targetTileX = (int) targetX/ Utilities.TILE_SIZE;
        targetTileY = (int) targetY / Utilities.TILE_SIZE;
    
        dangerMap = new boolean[width][height];
        for (Projectile p : projectiles) {
            if(p.getOwner() != this) {
                double currpx = p.getX();
                double currpy = p.getY();
                double px = p.getX();
                double py = p.getY();
                double angle = p.getAngle();
                double dx = Math.cos(angle);
                double dy = Math.sin(angle);

                for (int t = 0; t < 20; t++) {
                    int tx = (int)(px / Utilities.TILE_SIZE);
                    int ty = (int)(py / Utilities.TILE_SIZE);

                    double tileCenterX = tx * Utilities.TILE_SIZE + Utilities.TILE_SIZE / 2.0;
                    double tileCenterY = ty * Utilities.TILE_SIZE + Utilities.TILE_SIZE / 2.0;
                    double myDist = (Math.abs(getX() + (lastX == -1 || lastY == -1 ? Utilities.TILE_SIZE : 0) - tileCenterX) + Math.abs(getY() + (lastX == -1 || lastY == -1 ? Utilities.TILE_SIZE : 0) - tileCenterY))/Utilities.TILE_SIZE;
                    int mySteps = (int)Math.ceil(myDist / (getSpeed()));
                    double projDist = Math.sqrt((currpx - tileCenterX) * (currpx - tileCenterX) + (currpy - tileCenterY) * (currpy - tileCenterY))/Utilities.TILE_SIZE;
                    int projSteps = (int)Math.ceil(projDist / (p.getProjectileSpeed()));

                    if (tx >= 0 && tx < width && ty >= 0 && ty < height) {
                        if(tiles[ty][tx] == Utilities.WALL) {
                            break;
                        }
                        if (mySteps > projSteps - (p.getProjectileSpeed() - getSpeed())/3 && mySteps < projSteps + (p.getProjectileSpeed() - getSpeed())/3) {
                            dangerMap[tx][ty] = true;
                        }
                    }
                    px += dx * Utilities.TILE_SIZE;
                    py += dy * Utilities.TILE_SIZE;
                }
            }
        }

        for(int i = 0; i < width; i++) {
            for(int j = 0; j < height; j++) {
                if(tiles[j][i] == Utilities.WALL || dangerMap[i][j]) {  
                    visited[i][j] = true;
                }else{
                    visited[i][j] = false;
                }
                if(lastX == -1 || lastY == -1) {
                    if(robotTileXPlus == i && robotTileYPlus == j) {
                        distances[i][j] = 0;
                    }else{
                        distances[i][j] = Integer.MAX_VALUE;
                    }
                }else {
                    if(robotTileX == i && robotTileY == j) {
                        distances[i][j] = 0;
                    }else{
                        distances[i][j] = Integer.MAX_VALUE;
                    }
                }
            }
        }

        for(Robot robot : otherRobots) {
            for(int i = Math.max(robot.getX()/Utilities.TILE_SIZE - distancesDist + 1, 0); i < Math.min(robot.getX()/Utilities.TILE_SIZE + distancesDist, tiles[0].length); i++) {
                for(int j = Math.max(robot.getY()/Utilities.TILE_SIZE - distancesDist + 1, 0); j < Math.min(robot.getY()/Utilities.TILE_SIZE + distancesDist, tiles.length); j++) {
                    visited[i][j] = true;
                }
            }
        }

        if(dangerMap[(int) targetX/Utilities.TILE_SIZE][(int) targetY/Utilities.TILE_SIZE]) {
            double[] newTarget = newTarget(targetX, targetY);
            targetX = newTarget[0];
            targetY = newTarget[1];
        }

        Queue<Tile> priority = new LinkedList<>();
        if(lastX == -1 || lastY == -1) {
            priority.add(new Tile(robotTileXPlus, robotTileYPlus));
            visited[robotTileXPlus][robotTileYPlus] = true;
        } else {
            priority.add(new Tile(robotTileX, robotTileY));
            visited[robotTileX][robotTileY] = true;
        }
        while(!priority.isEmpty()) {
            int curr = priority.size();
            for(int i = 0; i < curr; i++) {
                Tile tile = priority.poll();
                int tileX = tile.tileX();
                int tileY = tile.tileY();
                if(tileX + 1 < visited.length && !visited[tileX + 1][tileY]) {
                    visited[tileX + 1][tileY] = true;
                    distances[tileX + 1][tileY] = distances[tileX][tileY] + (tiles[tileY][tileX + 1] == Utilities.MUD ? 2 : 1);
                    priority.add(new Tile(tileX + 1, tileY));
                }
                if(tileX - 1 >= 0 && !visited[tileX - 1][tileY]) {
                    visited[tileX - 1][tileY] = true;
                    distances[tileX - 1][tileY] = distances[tileX][tileY] + (tiles[tileY][tileX - 1] == Utilities.MUD ? 2 : 1);
                    priority.add(new Tile(tileX - 1, tileY));
                }
                if(tileY + 1 < visited[0].length && !visited[tileX][tileY + 1]) {
                    visited[tileX][tileY + 1] = true;
                    distances[tileX][tileY + 1] = distances[tileX][tileY] + (tiles[tileY + 1][tileX] == Utilities.MUD ? 2 : 1);
                    priority.add(new Tile(tileX, tileY + 1));
                }
                if(tileY - 1 >= 0 && !visited[tileX][tileY - 1]) {
                    visited[tileX][tileY - 1] = true;
                    distances[tileX][tileY - 1] = distances[tileX][tileY] + (tiles[tileY - 1][tileX] == Utilities.MUD ? 2 : 1);
                    priority.add(new Tile(tileX, tileY - 1));
                }
            }
        }

        int x = (int) targetTileX;
        int y = (int) targetTileY;
        int currDist = distances[x][y];
        while(currDist != 1) {
            int min = Math.min(Math.min(x + 1 < distances.length ? distances[x + 1][y] : Integer.MAX_VALUE, x - 1 >= 0 ? distances[x - 1][y] : Integer.MAX_VALUE), Math.min(y + 1 < distances[0].length ? distances[x][y + 1] : Integer.MAX_VALUE, y - 1 >= 0 ? distances[x][y - 1] : Integer.MAX_VALUE));
            if(x + 1 < distances.length && distances[x + 1][y] == min) {
                x += 1;
            }else if(x - 1 >= 0 && distances[x - 1][y] == min) {
                x -= 1;
            }else if(y + 1 < distances[0].length && distances[x][y + 1] == min) {
                y += 1;
            }else if(y - 1 >= 0 && distances[x][y - 1] == min) {
                y -= 1;
            }
            currDist = distances[x][y];
            if(currDist == Integer.MAX_VALUE) {
                runAway(projectiles);
                return;
            }
        }
        if(lastX == -1 || lastY == -1) {
            xMovement = x - robotTileXPlus;
            yMovement = y - robotTileYPlus;
            lastX = xMovement;
            lastY = yMovement;
        } else {
            xMovement = x - robotTileX;
            yMovement = y - robotTileY;
            lastX = x - robotTileX;
            lastY = y - robotTileY;
        }
    } 
    
    public void takeDamage(int amount) {
        System.out.println("Took Damage");
        super.takeDamage(amount);
    }

    public boolean[][] getDangerMap() {
        return dangerMap;
    }

    private double[] newTarget(double targetX, double targetY) {
        double[] targets = new double[2];
        boolean[][] already = new boolean[distances.length][distances[0].length];
        targets[0] = targetX;
        targets[1] = targetY;

        Queue<double[]> priority = new LinkedList<>();
        priority.add(targets);
        already[(int) targetX/Utilities.TILE_SIZE][(int) targetY/Utilities.TILE_SIZE] = true;
        while(!priority.isEmpty()) {
            double[] target = priority.poll();
            if(distances[Math.min((int) target[0]/Utilities.TILE_SIZE + 1, distances.length - 1)][(int) target[1]/Utilities.TILE_SIZE] != Integer.MAX_VALUE) {
                double[] next = new double[2];
                next[0] = target[0] + Utilities.TILE_SIZE;
                next[1] = target[1];
                return next;
            }else if(!already[Math.min((int) target[0]/Utilities.TILE_SIZE + 1, distances.length - 1)][(int) target[1]/Utilities.TILE_SIZE]) {
                double[] next = new double[2];
                next[0] = target[0] + Utilities.TILE_SIZE;
                next[1] = target[1];
                priority.add(next);
                already[Math.min((int) target[0]/Utilities.TILE_SIZE + 1, distances.length - 1)][(int) target[1]/Utilities.TILE_SIZE] = true;
            }
            if(distances[Math.max((int) target[0]/Utilities.TILE_SIZE - 1, 0)][(int) target[1]/Utilities.TILE_SIZE] != Integer.MAX_VALUE) {
                double[] next = new double[2];
                next[0] = target[0] - Utilities.TILE_SIZE;
                next[1] = target[1];
                return next;
            }else if(!already[Math.max((int) target[0]/Utilities.TILE_SIZE - 1, 0)][(int) target[1]/Utilities.TILE_SIZE]) {
                double[] next = new double[2];
                next[0] = target[0] - Utilities.TILE_SIZE;
                next[1] = target[1];
                priority.add(next);
                already[Math.max((int) target[0]/Utilities.TILE_SIZE - 1, 0)][(int) target[1]/Utilities.TILE_SIZE] = true;
            }
            if(distances[(int) target[0]/Utilities.TILE_SIZE][Math.min((int) target[1]/Utilities.TILE_SIZE + 1, distances[0].length - 1)] != Integer.MAX_VALUE) {
                double[] next = new double[2];
                next[1] = target[1] + Utilities.TILE_SIZE;
                next[0] = target[0];
                return next;
            }else if(!already[(int) target[0]/Utilities.TILE_SIZE][Math.min((int) target[1]/Utilities.TILE_SIZE + 1, distances[0].length - 1)]) {
                double[] next = new double[2];
                next[0] = target[0];
                next[1] = target[1] + Utilities.TILE_SIZE;
                priority.add(next);
                already[(int) target[0]/Utilities.TILE_SIZE][Math.min((int) target[1]/Utilities.TILE_SIZE + 1, distances[0].length - 1)] = true;
            }
            if(distances[(int) target[0]/Utilities.TILE_SIZE][Math.max((int) target[1]/Utilities.TILE_SIZE - 1, 0)] != Integer.MAX_VALUE) {
                double[] next = new double[2];
                next[1] = target[1] - Utilities.TILE_SIZE;
                next[0] = target[0];
                return next;
            }else if(!already[(int) target[0]/Utilities.TILE_SIZE][Math.max((int) target[1]/Utilities.TILE_SIZE - 1, 0)]) {
                double[] next = new double[2];
                next[0] = target[0];
                next[1] = target[1] - Utilities.TILE_SIZE;
                priority.add(next);
                already[(int) target[0]/Utilities.TILE_SIZE][Math.max((int) target[1]/Utilities.TILE_SIZE - 1, 0)] = true;
            }
        }
        return new double[2];
    }

    private void runAway(ArrayList<Projectile> projectiles) {
        double[] moveScores = new double[4]; // up, down, left, right
        int currX = getX();
        int currY = getY();

        for (Projectile p : projectiles) {
            double px = p.getX();
            double py = p.getY();
            double distUp = Math.hypot(currX - px, (currY - Utilities.TILE_SIZE) - py);
            moveScores[0] -= 1.0 / (distUp + 1);
            moveScores[2] += 1.0 / ((distUp + 1) * 2);
            moveScores[3] += 1.0 / ((distUp + 1) * 2);
            double distDown = Math.hypot(currX - px, (currY + Utilities.TILE_SIZE) - py);
            moveScores[1] -= 1.0 / (distDown + 1);
            moveScores[2] += 1.0 / ((distDown + 1) * 2);
            moveScores[3] += 1.0 / ((distDown + 1) * 2);
            double distLeft = Math.hypot((currX - Utilities.TILE_SIZE) - px, currY - py);
            moveScores[2] -= 1.0 / (distLeft + 1);
            moveScores[0] += 1.0 / ((distDown + 1) * 2);
            moveScores[1] += 1.0 / ((distDown + 1) * 2);
            double distRight = Math.hypot((currX + Utilities.TILE_SIZE) - px, currY - py);
            moveScores[3] -= 1.0 / (distRight + 1);
            moveScores[0] += 1.0 / ((distDown + 1) * 2);
            moveScores[1] += 1.0 / ((distDown + 1) * 2);
        }

        int bestDir = 0;
        for (int i = 1; i < 4; i++) {
            if (moveScores[i] > moveScores[bestDir]) {
                if(distances[getX()/Utilities.TILE_SIZE][getY()/Utilities.TILE_SIZE + 1] == Utilities.WALL && i == 0) {
                    continue;
                }
                if(distances[getX()/Utilities.TILE_SIZE][getY()/Utilities.TILE_SIZE - 1] == Utilities.WALL && i == 1) {
                    continue;
                }
                if(distances[getX()/Utilities.TILE_SIZE + 1][getY()/Utilities.TILE_SIZE] == Utilities.WALL && i == 3) {
                    continue;
                }
                if(distances[getX()/Utilities.TILE_SIZE - 1][getY()/Utilities.TILE_SIZE] == Utilities.WALL && i == 2) {
                    continue;
                }
                bestDir = i;
            }
        }
        switch (bestDir) {
            case 0: 
                xMovement = 0;
                yMovement = -1;
                break;
            case 1: 
                xMovement = 0;
                yMovement = 1;
                break;
            case 2:
                xMovement = -1;
                yMovement = 0;
                break;
            case 3:
                xMovement = 1;
                yMovement = 0;
                break;
        }
    }

    private boolean hasLineOfSight(Robot target) {
        int startX = getX() + Utilities.ROBOT_SIZE / 2;
        int startY = getY() + Utilities.ROBOT_SIZE / 2;
        int targetLeft = target.getX();
        int targetRight = target.getX() + Utilities.ROBOT_SIZE - 1;
        int targetTop = target.getY();
        int targetBottom = target.getY() + Utilities.ROBOT_SIZE - 1;

        int[][] tiles = this.tiles;
        for (int tx = targetLeft; tx <= targetRight; tx += Math.max(1, Utilities.ROBOT_SIZE / 2)) {
            for (int ty = targetTop; ty <= targetBottom; ty += Math.max(1, Utilities.ROBOT_SIZE / 2)) {
                if (lineClear(startX, startY, tx, ty, tiles)) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean lineClear(int x0, int y0, int x1, int y1, int[][] tiles) {
        int dx = Math.abs(x1 - x0);
        int dy = Math.abs(y1 - y0);
        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;
        int err = dx - dy;
        while (true) {
            int tileX = x0 / Utilities.TILE_SIZE;
            int tileY = y0 / Utilities.TILE_SIZE;
            if (tiles[tileY][tileX] == Utilities.WALL) return false;
            if (x0 == x1 && y0 == y1) break;
            int e2 = 2 * err;
            if (e2 > -dy) { err -= dy; x0 += sx; }
            if (e2 < dx) { err += dx; y0 += sy; }
        }
        return true;
    }
}