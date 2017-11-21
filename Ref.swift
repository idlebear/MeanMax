//
//  Ref.swift
//  WaterPigs
//
//  Ported from the competition referee at: https://github.com/CodinGame/MeanMax/blob/master/Referee.java
//  Created by Barry Gilhuly on 2017-11-19.
//  Copyright © 2017 Red Maples Consulting. All rights reserved.
//

#if os(Linux)
    import Glibc
#endif
import Foundation

public struct StderrOutputStream: TextOutputStream {
    public mutating func write(_ string: String) { fputs(string, stderr) }
}
public var errStream = StderrOutputStream()

/**
 
 Point -- a cartesian representation of a point on a plane - (x,y)
 
 */
public
struct Point : CustomStringConvertible, Codable {
    
    var _x: Double
    var _y: Double
    
    public init( _ x: Int, _ y: Int ) {
        self.init( Double(x), Double(y))
    }
    
    public init( _ x: CGFloat, _ y: CGFloat ) {
        self.init( Double(x), Double(y))
    }
    
    public init( _ x: Double, _ y: Double ) {
        self._x = x
        self._y = y
    }
    
    public var x:Double {
        get {
            return self._x
        }
        set {
            self._x = newValue
        }
    }
    
    public var y:Double {
        get {
            return self._y
        }
        set {
            self._y = newValue
        }
    }
    
    public func distance(to point: Point) -> Double {
        let dx = self.x - point.x
        let dy = self.y - point.y
        return sqrt(dx * dx + dy * dy)
    }
    
    public func walkingDistance(to point: Point) -> Double {
        let dx = self.x - point.x
        let dy = self.y - point.y
        return abs(dx) + abs(dy)
    }
    
    public var description: String {
        return "(\(x),\(y))"
    }
    
}

extension Point: Equatable {
    public static func == ( lhs: Point, rhs: Point ) -> Bool {
        return ( lhs.x == rhs.x ) && ( lhs.y == rhs.y )
    }
}


extension Point {
    public static func centre(of points:[Point]) -> Point {
        guard points.count > 0 else { return Point(0,0) }
        
        var midX = 0.0
        var midY = 0.0
        for point in points {
            midX += point.x
            midY += point.y
        }
        midX /= Double(points.count)
        midY /= Double(points.count)
        return Point(midX, midY)
    }
    
    public static func dominantAxis(of points: [Point]) -> (Point,Double) {
        let centre = Point.centre(of: points)
        
        // Calculate the perpendicular (as opposed to the vertical) minimized difference
        // to the line of best fit -- based on a paper found at:
        //
        // https://www.engr.colostate.edu/~dga/dga/papers/least_squares.pdf
        
        var bPrime = 0.0
        var sumX = 0.0
        var sumY = 0.0
        for point in points {
            let xDiff = point.x - centre.x
            let yDiff = point.y - centre.y
            sumX += xDiff * xDiff
            sumY += yDiff * yDiff
            bPrime += xDiff * yDiff
        }
        let aPrime = sumX - sumY
        
        if bPrime == 0 {
            if sumX == 0 {
                return (centre, 0) // horizontal line
            }
            return (centre, Double.pi/2)  // vertical line
        } else {
            let A = Double(2 * bPrime)
            let B = -1 * ( Double(aPrime) + sqrt( Double(aPrime * aPrime + 4 * bPrime * bPrime) ) )
            // don't need the y intercept
            // let CC = AA * Double(centre.x) + BB * Double(centre.y)
            
            // equation for line A * x + B * y = C, so slope == -A/B
            return (centre, atan(-1*Double(A)/Double(B)))
        }
    }
}

extension Point: Hashable {
    // Codingame specific functions
    
    // Move the point to an other point for a given distance
    func move( towards p: Point, for dist: Double ) -> Point {
        let diff = distance(to: p)
        if diff < EPSILON {
            return self
        }
        
        let dx = p.x - x
        let dy = p.y - y
        let coef = dist / diff
        
        return Point( x + dx * coef, y + dy * coef )
    }
    
    func isInRange( _ p: Point, _ range: Double) -> Bool {
        return (p.x != x && p.y != y) && distance(to:p) <= range
    }
    
    public
    var hashValue: Int {
        let prime = 31
        var result = 1
        let tx = x.bitPattern
        result = prime * result + Int(tx ^ (tx>>32))
        let ty = y.bitPattern
        result = prime * result + Int(ty ^ (ty>>32))
        return result
    }
    
}

public struct Polar {
    var angle: Double
    var length: Double
    
    init( _ angle: Double, _ length: Double ) {
        self.angle = angle
        self.length = length
    }
    
    init( _ point: Point ) {
        angle = atan2( Double(point.y), Double(point.x) )
        if angle < 0 {
            angle += Double.pi
        }
        length = sqrt( Double(point.x*point.x + point.y*point.y))
    }
    
    func toPoint() -> Point {
        let x = length * cos(angle)
        let y = length * sin(angle)
        return Point(x,y)
    }
    
    func rotate( by angle: Double) -> Polar {
        return Polar( self.angle + angle, self.length )
    }
}

// MARK: Vector Implementation
// *****************************************************************************


public
struct Vector: Codable {
    
    public var dx: Double
    public var dy: Double
    
    public init( dx: Double, dy: Double ) {
        self.dx = dx
        self.dy = dy
    }
    
    public init( dx: Int, dy: Int ) {
        self.init( dx: Double(dx), dy: Double(dy))
    }
    
    public init( dx: CGFloat, dy: CGFloat ) {
        self.init( dx: Double(dx), dy: Double(dy))
    }
    
    public init( start: Point, end: Point ) {
        self.init( dx: end.x - start.x, dy: end.y - start.y )
    }
    
    public init(length: Double, angle: Double ) {
        let dx = length * cos(angle)
        let dy = length * sin(angle)
        self.init( dx: dx, dy: dy )
    }
    
    public init(length: CGFloat, angle: CGFloat) {
        self.init( length: Double(length), angle: Double(angle) )
    }
    
    public var length: Double {
        return sqrt( length_squared )
    }
    
    public var length_squared: Double {
        return dx * dx + dy * dy
    }
    
    public func asUnitVector() -> Vector {
        guard length > 0 else { return Vector(dx: 0,dy: 0) }
        
        return Vector( dx: (dx/length), dy: (dy/length) )
    }
    
    public func getRightNormal() -> Vector {
        return Vector( dx: dy, dy: -dx ).asUnitVector()
    }
    
    public func getLeftNormal() -> Vector {
        return Vector( dx: -dy, dy: dx ).asUnitVector()
    }
    
    public func reversed() -> Vector {
        return Vector( dx: -dx, dy: -dy )
    }
    
    public func flippedHorizontal() -> Vector {
        return Vector(dx: -dx, dy: dy)
    }
    
    public func flippedVertical() -> Vector {
        return Vector(dx: dx, dy: -dy)
    }
    
    public func slope() -> Double? {
        guard dx != 0 else { return nil }
        return dy / dx
    }
    
    public func scaled(by scaleFactor: Double) -> Vector {
        
        //        cos( theta ) = x / length
        //        sin( theta ) = y / length
        //
        //        new length = length * scaleFactor = new_x / cos( theta )
        //                                          = new_x / ( x / length )
        //                                          = ( new_x * length ) / x
        //
        //        therefore,   new x = ( scaleFactor * x * length ) / length
        //                           =   scaleFactor * x
        
        return Vector(dx: dx*scaleFactor, dy: dy*scaleFactor)
    }
    
    public func angle() -> Double {
        return atan2(dy, dx)
    }
    
}

// DOT product
infix operator !
public func ! (lhs: Vector, rhs: Vector) -> Double {
    return lhs.dx * rhs.dx + lhs.dy * rhs.dy
}

public func == (lhs: Vector, rhs: Vector) -> Bool {
    if lhs.dx == rhs.dx && lhs.dy == rhs.dy {
        return true
    }
    return false
}

// CROSS product -- the cross product of two 2d vectors doesn't exist
// in 2D space in the same way as in 3D (where it results in a vector perpendicular
// to the two input vectors).  However, if one assumes that the third term (Z),
// is zero, then the cross product resolves to the Determinant of the two
// vectors:
//           x1 * y2 - y1 * x2
// and happens to be the area of the parallelogram defined by the input
// vectors
public func * (lhs: Vector, rhs: Vector) -> Double {
    return lhs.dx * rhs.dy - lhs.dy * rhs.dx
}

public func + (lhs: Vector, rhs: Vector) -> Vector {
    return Vector(dx: lhs.dx + rhs.dx, dy: lhs.dy + rhs.dy)
}

public func - (lhs: Vector, rhs: Vector) -> Vector {
    return Vector(dx: lhs.dx - rhs.dx, dy: lhs.dy - rhs.dy)
}

public func * (lhs: Vector, rhs: Double ) -> Vector {
    return lhs.scaled(by: rhs)
}

public func * (lhs: Double, rhs: Vector) -> Vector {
    return rhs.scaled(by: lhs)
}

public func + (lhs: Point, rhs: Vector) -> Point {
    return Point(lhs.x + rhs.dx, lhs.y + rhs.dy)
}


// Protocol and extension to get the sign from a signed numerical value
protocol Signable {
    init()
    init(_:Int)
    static func <(lhs:Self,rhs:Self) -> Bool
}

extension Signable where Self:SignedNumeric {
    func sign() -> Self {
        return (self < Self()) ? Self(-1) : Self(1)
    }
}

extension Int: Signable {}
extension Int32: Signable {}
extension Int64: Signable {}
extension Double: Signable {}

protocol Numeric {
    var asDouble: Double { get }
    init(_: Double)
}

extension Int: Numeric {var asDouble: Double {return Double(self)}}
extension Int64: Numeric {var asDouble: Double {return Double(self)}}
extension Float: Numeric {var asDouble: Double {return Double(self)}}
extension Double: Numeric {var asDouble: Double {return Double(self)}}
extension CGFloat: Numeric {var asDouble: Double {return Double(self)}}

extension Array where Element: Numeric {
    var mean: Element {
        get {
            let mean = self.reduce(0.0){($0 + $1.asDouble)}
            let n = Double(self.count)
            return Element(mean/n)
        }
    }
    
    var standardDeviation : Element {
        get {
            let sss = self.reduce((0.0, 0.0)){ return ($0.0 + $1.asDouble, $0.1 + ($1.asDouble * $1.asDouble))}
            let n = Double(self.count)
            return Element(sqrt(sss.1/n - (sss.0/n * sss.0/n)))
        }
    }
}

/****** Timer ********
 */

public struct Timer {
    var startData:timespec = timespec()
    var endData:timespec = timespec()
    
    public init() {
        start()
    }
    
    public mutating func start() {
        clock_gettime( CLOCK_MONOTONIC, &startData )
        endData = startData
    }
    
    public mutating func end() -> Int64 {
        clock_gettime( CLOCK_MONOTONIC, &endData )
        return result()
    }
    
    public func result() -> Int64 {
        return Int64(endData.tv_sec - startData.tv_sec) * 1_000_000_000 + Int64(endData.tv_nsec - startData.tv_nsec)
    }
    
    public func peek() -> Int64 {
        var tempData: timespec = timespec()
        clock_gettime( CLOCK_MONOTONIC, &tempData )
        return Int64(tempData.tv_sec - startData.tv_sec) * 1_000_000_000 + Int64(tempData.tv_nsec - startData.tv_nsec)
    }
    
    public func peek_uSec() -> Double {
        return Double(peek())/1_000.0
    }
    
    public func peek_mSec() -> Double {
        return Double(peek()/1_000)/1_000.0
    }
    
    public static func benchmark( _ fn:() -> Void ) -> (Double,Double) {
        var recordedTimes:[Double] = []
        for _ in 0..<100 {
            let t:Timer = Timer()
            fn()
            recordedTimes.append(t.peek_mSec())
            if recordedTimes.count % 10 == 0 {
                let mean = recordedTimes.mean
                let standardDeviation = recordedTimes.standardDeviation
                print( "\(recordedTimes.count): \(mean) \(standardDeviation)" )
            }
        }
        let mean = recordedTimes.mean
        let standardDeviation = recordedTimes.standardDeviation
        return( mean, standardDeviation )
    }
}

///
///  Xoroshiro128+ implementation ported from http://xoroshiro.di.unimi.it/#shootout
///  Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
///
///  Splitmx64 implementation ported from http://xoroshiro.di.unimi.it/splitmix64.c
///  Written in 2015 by Sebastiano Vigna (vigna@acm.org)
///
///  Adapted components from http://cocoawithlove.com to implement simple (primarily smaller)
///  random number generator.  Created by Matt Gallagher on 2016/05/17.
///
///

final
public class Random {
    
    typealias StateType = (UInt64,UInt64)
    
    static var seed: UInt64 = 0
    static var state: StateType? = nil
    
    class FileDescriptor {
        let value: CInt
        init() {
            value = open("/dev/urandom", O_RDONLY)
            precondition(value >= 0)
        }
        deinit {
            close(value)
        }
    }
    
    
    public init( seed:UInt64 ) {
        if Random.state == nil {
            Random.seed = seed
            Random.state = (Random.nextSeed(),Random.nextSeed())
        }
    }
    
    public convenience init() {
        var seed: UInt64 = 0
        if Random.state == nil {
            let fd: FileDescriptor = FileDescriptor()
            _ = withUnsafeMutablePointer(to: &seed, {
                $0.withMemoryRebound(to: UInt8.self, capacity: MemoryLayout<UInt64>.size, {
                    return read(fd.value, $0, MemoryLayout<UInt64>.size)
                })
            })
        }
        self.init(seed: seed)
    }
    
    static func seed(with seed: UInt64) -> Void {
        Random.seed = seed
        Random.state = (nextSeed(),nextSeed())
    }
    
    fileprivate
    static func nextSeed() -> UInt64 {
        seed = seed &+ (0x9E3779B97F4A7C15 as UInt64)
        var z: UInt64 = seed
        z = (z ^ (z >> 30)) &* (0xBF58476D1CE4E5B9 as UInt64)
        z = (z ^ (z >> 27)) &* (0x94D049BB133111EB as UInt64)
        return z ^ (z >> 31)
    }
    
    fileprivate static func rotl( _ x: UInt64, _ k: UInt64 ) -> UInt64 {
        return (x << k) | (x >> (64 - k))
    }
    
    fileprivate
    static func next() -> UInt64 {
        guard Random.state != nil else {
            return 0
        }
        let s0 = Random.state!.0
        var s1 = Random.state!.1
        
        let result = s0 &+ s1
        s1 ^= s0
        Random.state!.0 = Random.rotl(s0, 55) ^ s1 ^ (s1<<14) // a, b
        Random.state!.1 = Random.rotl(s1, 36)                 // c
        
        return result
    }
    
    
    final
    public func double() -> Double {
        return Double(Random.next() & 0x001f_ffff_ffff_ffff) * (1.0 / 9007199254740992.0)
    }
    
    final
    public func int( from low: Int = 0, to high: Int ) -> Int {
        guard low < high else { return 0 }
        
        let diff = UInt64( high - low )
        let limit = UInt64.max & diff
        
        var result: UInt64
        repeat {
            result = Random.next()
        } while result < limit
        result = result % diff
        
        return low + Int(result)
    }
}


let GAME_VERSION = 3;

let SPAWN_WRECK = false;
let LOOTER_COUNT = 3;
let REAPER_SKILL_ACTIVE = true;
let DESTROYER_SKILL_ACTIVE = true;
let DOOF_SKILL_ACTIVE = true;

let MAP_RADIUS = 6000.0;
var TANKERS_BY_PLAYER_MIN = 1;
var TANKERS_BY_PLAYER_MAX = 3;

let WATERTOWN_RADIUS = 3000.0;

let TANKER_THRUST = 500;
let TANKER_EMPTY_MASS = 2.5;
let TANKER_MASS_BY_WATER = 0.5;
let TANKER_FRICTION = 0.40;
let TANKER_RADIUS_BASE = 400.0;
let TANKER_RADIUS_BY_SIZE = 50.0;
let TANKER_EMPTY_WATER = 1;
let TANKER_MIN_SIZE = 4;
let TANKER_MAX_SIZE = 10;
let TANKER_MIN_RADIUS = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * Double(TANKER_MIN_SIZE);
let TANKER_MAX_RADIUS = TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * Double(TANKER_MAX_SIZE);
let TANKER_SPAWN_RADIUS = 8000.0;
let TANKER_START_THRUST = 2000;

let MAX_THRUST = 300;
let MAX_RAGE = 300;
let WIN_SCORE = 50;

let REAPER_MASS = 0.5;
let REAPER_FRICTION = 0.20;
let REAPER_SKILL_DURATION = 3;
let REAPER_SKILL_COST = 30;
let REAPER_SKILL_ORDER = 0;
let REAPER_SKILL_RANGE = 2000.0;
let REAPER_SKILL_RADIUS = 1000.0;
let REAPER_SKILL_MASS_BONUS = 10.0;

let DESTROYER_MASS = 1.5;
let DESTROYER_FRICTION = 0.30;
let DESTROYER_SKILL_DURATION = 1;
let DESTROYER_SKILL_COST = 60;
let DESTROYER_SKILL_ORDER = 2;
let DESTROYER_SKILL_RANGE = 2000.0;
let DESTROYER_SKILL_RADIUS = 1000.0;
let DESTROYER_NITRO_GRENADE_POWER = 1000;

let DOOF_MASS = 1.0;
let DOOF_FRICTION = 0.25;
let DOOF_RAGE_COEF = 1.0 / 100.0;
let DOOF_SKILL_DURATION = 3;
let DOOF_SKILL_COST = 30;
let DOOF_SKILL_ORDER = 1;
let DOOF_SKILL_RANGE = 2000.0;
let DOOF_SKILL_RADIUS = 1000.0;

let LOOTER_RADIUS = 400.0;
let LOOTER_REAPER = 0;
let LOOTER_DESTROYER = 1;
let LOOTER_DOOF = 2;

let TYPE_TANKER = 3;
let TYPE_WRECK = 4;
let TYPE_REAPER_SKILL_EFFECT = 5;
let TYPE_DOOF_SKILL_EFFECT = 6;
let TYPE_DESTROYER_SKILL_EFFECT = 7;

let EPSILON = 0.00001;
let MIN_IMPULSE = 30.0;
let IMPULSE_COEFF = 0.5;

// Center of the map
var WATERTOWN = Point(0, 0);

protocol IdGenerator: IteratorProtocol where Element == Int {
    mutating func next() -> Self.Element?
}

struct GlobalId: IdGenerator {
    typealias Element = Int
    var id: Int = 0
    public mutating func next() -> Int? {
        defer { id += 1 }
        return id
    }
}


enum Move {
    case wait, accel( Point, Int ), skill( Point )
    
    var description: String {
        switch self {
        case .wait:
            return "WAIT"
        case .accel( let point, let power ):
            return "\(point.x) \(point.y) \(power)"
        case .skill( let target ):
            return "SKILL \(target.x) \(target.y)"
        }
    }
}


enum UnitType {
    case reaper, destroyer, doof, tanker, wreck, tar, oil, unknown
    
    init( _ val: Int ) {
        switch  val {
        case 0:
            self = .reaper
        case 1:
            self = .destroyer
        case 2:
            self = .doof
        case 3:
            self = .tanker
        case 4:
            self = .wreck
        case 5:
            self = .tar
        case 6:
            self = .oil
        default:
            self = .unknown
        }
    }
    
    func friction() -> Double {
        switch self {
        case .reaper:
            return 0.2
        case .destroyer:
            return 0.3
        case .doof:
            return 0.25
        case .tanker:
            return 0.2
        default:
            return 1
        }
    }
    
    var description: String {
        switch self {
        case .reaper:
            return "Reaper"
        case .destroyer:
            return "Destroyer"
        case .doof:
            return "Doof"
        case .tanker:
            return "Tanker"
        case .wreck:
            return "Wreck"
        case .tar:
            return "TarSpill"
        case .oil:
            return "OilSlick"
        case .unknown:
            return "Unknown"
        }
    }
    
    var shortName: String {
        switch self {
        case .reaper:
            return "R"
        case .destroyer:
            return "D"
        case .doof:
            return "d"
        case .tanker:
            return "T"
        case .wreck:
            return "W"
        case .tar:
            return "T"
        case .oil:
            return "O"
        case .unknown:
            return "U"
        }
    }
}



class Ref {
    
    var TANKERS_BY_PLAYER: Int

    // Global first free id for all elements on the map
    static var GLOBAL_ID = GlobalId()
    
    // The null collision
    static var NULL_COLLISION:Collision = Collision(t:1.0 + EPSILON);

    class Unit : Hashable {
        
        let id: Int
        let type: Int
        var pos: Point
        let radius: Double
        var vx: Double = 0
        var vy: Double = 0
        var mass: Double = 0
        var friction: Double = 0
        var known: Bool = false
        
        init( type: Int, pos: Point, radius: Double ) {
            self.id = Ref.GLOBAL_ID.next()!
            self.type = type
            self.pos = pos
            self.radius = radius
        }
        

        func move(_ t: Double) {
            pos.x += vx * t;
            pos.y += vy * t;
        }
        
        func speed() -> Double {
            return sqrt(vx * vx + vy * vy)
        }
        
        final
        var hashValue: Int {
            let prime = 31;
            var result = 1;
            result = prime * result + id;
            return result;
        }
        
        static func ==(lhs: Ref.Unit, rhs: Ref.Unit) -> Bool {
            if lhs === rhs {
                return true
            }
            if lhs.id != rhs.id {
                return false
            }
            return true
        }
        
        func getFrameId() -> String {
            return "\(id)"
        }
        
        func toFrameData() -> String {
            if (known) {
                return "\(getFrameId()) \(round(pos.x)) \(round(pos.y)) \(round(vx)) \(round(vy))"
            }
            known = true
            return "\(getFrameId()) \(round(pos.x)) \(round(pos.y)) \(round(vx)) \(round(vy)) \(type) \(round(radius))"
        }
        
        func thrust( _ p: Point, _ power: Int ) {
            let distance = pos.distance(to:p)
            
            // Avoid a division by zero
            if abs(distance) <= EPSILON {
                return
            }
            
            let coef = (Double(power) / mass) / distance
            vx += (p.x - pos.x) * coef;
            vy += (p.y - pos.y) * coef;
        }
        
        func isInDoofSkill(skillEffects: Set<SkillEffect>) -> Bool {
            for se in skillEffects {
                if se is DoofSkillEffect && pos.isInRange(se.pos, se.radius + radius) {
                    return true
                }
            }
            return false
        }
        
        func adjust(skillEffects: Set<SkillEffect> ) {
            pos.x = round(pos.x);
            pos.y = round(pos.y);
            
            if (isInDoofSkill(skillEffects: skillEffects)) {
                // No friction if we are in a doof skill effect
                vx = round(vx);
                vy = round(vy);
            } else {
                vx = round(vx * (1.0 - friction));
                vy = round(vy * (1.0 - friction));
            }
        }
        
        // Search the next collision with the map border
        func getCollision() -> Collision {
            // Check instant collision
            if pos.distance(to:WATERTOWN) + radius >= MAP_RADIUS {
                return Collision(t:0.0, a:self)
            }
            
            // We are not moving, we can't reach the map border
            if (vx == 0.0 && vy == 0.0) {
                return NULL_COLLISION
            }
            
            // Search collision with map border
            // Resolving: sqrt((x + t*vx)^2 + (y + t*vy)^2) = MAP_RADIUS - radius <=> t^2*(vx^2 + vy^2) + t*2*(x*vx + y*vy) + x^2 + y^2 - (MAP_RADIUS - radius)^2 = 0
            // at^2 + bt + c = 0;
            // a = vx^2 + vy^2
            // b = 2*(x*vx + y*vy)
            // c = x^2 + y^2 - (MAP_RADIUS - radius)^2
            
            let a = vx * vx + vy * vy;
            
            if a <= 0.0 {
                return NULL_COLLISION;
            }
            
            let b = 2.0 * (pos.x * vx + pos.y * vy);
            let c = pos.x * pos.x + pos.y * pos.y - (MAP_RADIUS - radius) * (MAP_RADIUS - radius);
            let delta = b * b - 4.0 * a * c;
            
            if (delta <= 0.0) {
                return NULL_COLLISION;
            }
            
            let t = (-b + sqrt(delta)) / (2.0 * a);
            
            if (t <= 0.0) {
                return NULL_COLLISION;
            }
            
            return Collision(t:t, a:self);
        }
        
        // Search the next collision with an other unit
        func getCollision(_ u: Unit) -> Collision {
            // Check instant collision
            if (pos.distance(to:u.pos) <= radius + u.radius) {
                return Collision(t:0.0, a:self, b:u);
            }
            
            // Both units are motionless
            if (vx == 0.0 && vy == 0.0 && u.vx == 0.0 && u.vy == 0.0) {
                return NULL_COLLISION;
            }
            
            // Change referencial
            // Unit u is not at point (0, 0) with a speed vector of (0, 0)
            let x2 = pos.x - u.pos.x;
            let y2 = pos.y - u.pos.y;
            let r2 = radius + u.radius;
            let vx2 = vx - u.vx;
            let vy2 = vy - u.vy;
            
            // Resolving: sqrt((x + t*vx)^2 + (y + t*vy)^2) = radius <=> t^2*(vx^2 + vy^2) + t*2*(x*vx + y*vy) + x^2 + y^2 - radius^2 = 0
            // at^2 + bt + c = 0;
            // a = vx^2 + vy^2
            // b = 2*(x*vx + y*vy)
            // c = x^2 + y^2 - radius^2
            
            let a = vx2 * vx2 + vy2 * vy2;
            
            if (a <= 0.0) {
                return NULL_COLLISION;
            }
            
            let b = 2.0 * (x2 * vx2 + y2 * vy2);
            let c = x2 * x2 + y2 * y2 - r2 * r2;
            let delta = b * b - 4.0 * a * c;
            
            if (delta < 0.0) {
                return NULL_COLLISION;
            }
            
            let t = (-b - sqrt(delta)) / (2.0 * a);
            
            if (t <= 0.0) {
                return NULL_COLLISION;
            }
            
            return Collision(t:t, a:self, b:u);
        }
        
        // Bounce between 2 units
        func bounce(_ u:Unit) {
            let mcoeff = (mass + u.mass) / (mass * u.mass);
            let nx = pos.x - u.pos.x;
            let ny = pos.y - u.pos.y;
            let nxnysquare = nx * nx + ny * ny;
            let dvx = vx - u.vx;
            let dvy = vy - u.vy;
            let product = (nx * dvx + ny * dvy) / (nxnysquare * mcoeff);
            var fx = nx * product;
            var fy = ny * product;
            let m1c = 1.0 / mass;
            let m2c = 1.0 / u.mass;
            
            vx -= fx * m1c;
            vy -= fy * m1c;
            u.vx += fx * m2c;
            u.vy += fy * m2c;
            
            fx = fx * IMPULSE_COEFF;
            fy = fy * IMPULSE_COEFF;
            
            // Normalize vector at min or max impulse
            let impulse = sqrt(fx * fx + fy * fy);
            var coeff = 1.0;
            if (impulse > EPSILON && impulse < MIN_IMPULSE) {
                coeff = MIN_IMPULSE / impulse;
            }
            
            fx = fx * coeff;
            fy = fy * coeff;
            
            vx -= fx * m1c;
            vy -= fy * m1c;
            u.vx += fx * m2c;
            u.vy += fy * m2c;
            
            let diff = (pos.distance(to:u.pos) - radius - u.radius) / 2.0;
            if (diff <= 0.0) {
                // Unit overlapping. Fix positions.
                pos = pos.move( towards:u.pos, for:diff - EPSILON)
                u.pos = u.pos.move( towards:pos, for:diff - EPSILON)
            }
        }
        
        // Bounce with the map border
        func bounce() {
            let mcoeff = 1.0 / mass;
            let nxnysquare = pos.x * pos.x + pos.y * pos.y;
            let product = (pos.x * vx + pos.y * vy) / (nxnysquare * mcoeff);
            var fx = pos.x * product;
            var fy = pos.y * product;
            
            vx -= fx * mcoeff;
            vy -= fy * mcoeff;
            
            fx = fx * IMPULSE_COEFF;
            fy = fy * IMPULSE_COEFF;
            
            // Normalize vector at min or max impulse
            let impulse = sqrt(fx * fx + fy * fy);
            var coeff = 1.0;
            if (impulse > EPSILON && impulse < MIN_IMPULSE) {
                coeff = MIN_IMPULSE / impulse;
            }
            
            fx = fx * coeff;
            fy = fy * coeff;
            vx -= fx * mcoeff;
            vy -= fy * mcoeff;
            
            let diff = pos.distance(to:WATERTOWN) + radius - MAP_RADIUS;
            if (diff >= 0.0) {
                // Unit still outside of the map, reposition it
                pos = pos.move(towards: WATERTOWN, for: diff + EPSILON)
            }
        }
        
        func getExtraInput() -> Int {
            return -1
        }
        
        func getExtraInput2() -> Int {
            return -1
        }
        
        func getPlayerIndex() -> Int {
            return -1
        }
    
    }
    
    class Wreck: Unit {
        var water : Int
        var player: Player? = nil
        
        init( x: Double, y: Double, water: Int, radius: Double ) {
            self.water = water;
            super.init( type: TYPE_WRECK, pos: Point(x,y), radius: radius )
        }
        
        override
        func  getFrameId() -> String {
            return "\(id)@\(water)"
        }
        
        override
        func toFrameData() -> String {
            if (known) {
                return getFrameId()
            }
            
            known = true;
            return "\(getFrameId()) \(round(pos.x)) \(round(pos.y)) 0 0 \(TYPE_WRECK) \(radius)"
        }
        
        // Reaper harvesting
        func harvest( players: [Player], skillEffects: Set<SkillEffect>) -> Bool {
            for p in players {
                let reaper = p.getReaper()!
                if pos.isInRange(reaper.pos, radius) && !reaper.isInDoofSkill(skillEffects: skillEffects) {
                    p.score += 1
                    water -= 1
                }
            }
            
            return water > 0
        }

    }
    
    class Tanker: Unit {
        var water: Int
        let size: Int
        
        var player:Player
        var killed:Bool = false
        
        init(size: Int, player:Player) {
            self.player = player
            self.size = size
            water = TANKER_EMPTY_WATER;
            super.init(type: TYPE_TANKER, pos: Point(0,0), radius: TANKER_RADIUS_BASE + TANKER_RADIUS_BY_SIZE * Double(size))
            
            mass = TANKER_EMPTY_MASS + TANKER_MASS_BY_WATER * Double(water)
            friction = TANKER_FRICTION;
        }
        
        override
        func getFrameId() -> String {
            return "\(id)@\(water)"
        }
        
        func die() -> Wreck? {
            // Don't spawn a wreck if our center is outside of the map
            if (pos.distance(to:WATERTOWN) >= MAP_RADIUS) {
                return nil
            }
            
            return  Wreck(x: pos.x, y:pos.y, water:water, radius:radius)
        }
        
        func isFull() -> Bool {
            return water >= size;
        }
        
        func play() {
            if (isFull()) {
                // Try to leave the map
                thrust(WATERTOWN, -TANKER_THRUST);
            } else if (pos.distance(to:WATERTOWN) > WATERTOWN_RADIUS) {
                // Try to reach watertown
                thrust(WATERTOWN, TANKER_THRUST);
            }
        }
        
        override func getCollision() -> Collision {
            // Tankers can go outside of the map
            return NULL_COLLISION;
        }
        
        override
        func getExtraInput() -> Int {
            return water;
        }
        
        override
        func getExtraInput2() -> Int {
            return size;
        }
    }
    
    enum LooterError : Error {
        case TooFarException, NoRageException
    }
    
    class Looter : Unit {
        let skillCost: Int
        let skillRange: Double
        var skillActive: Bool
        
        var player: Player
        
        var wantedThrustTarget: Point?
        var wantedThrustPower = 0
        var message: String?
        var attempt: Action?
        var skillResult: SkillResult?
        
        init(type: Int, player: Player, x: Double, y: Double, skillCost: Int, skillRange: Double, skillActive: Bool) {
            self.player = player;
            self.skillCost = skillCost
            self.skillRange = skillRange
            self.skillActive = skillActive
            super.init(type: type, pos: Point(x,y), radius: LOOTER_RADIUS)
        }
        
        // BUGBUG -- throws TooFarException, NoRageException
        func skill(_ p: Point) throws -> SkillEffect?  {
            guard player.rage >= skillCost else {
                throw LooterError.NoRageException
            }
            guard pos.distance(to: p) <= skillRange else {
                throw LooterError.TooFarException
            }
            player.rage -= skillCost;
            return skillImpl(p)
        }
        
        override
        func toFrameData() -> String {
            if (known) {
                return super.toFrameData();
            }
            return "\(super.toFrameData()) \(player.index)"
        }
        
        override
        func getPlayerIndex() -> Int {
            return player.index;
        }
        
        func skillImpl(_ p:Point) -> SkillEffect? {
            return nil
        }
        
        func setWantedThrust(_ target:Point, _ power:Int) {
            wantedThrustTarget = target;
            wantedThrustPower = max(0,min(power,MAX_THRUST))
        }
        
        func reset() {
            message = nil
            attempt = nil
            skillResult = nil
            wantedThrustTarget = nil
        }
        
    }
    
    class Reaper: Looter {
        init( player:Player, x:Double, y:Double ) {
            super.init(type: LOOTER_REAPER, player: player, x: x, y: y, skillCost: REAPER_SKILL_COST, skillRange: REAPER_SKILL_RANGE, skillActive: REAPER_SKILL_ACTIVE)
            mass = REAPER_MASS;
            friction = REAPER_FRICTION;
            skillActive = REAPER_SKILL_ACTIVE;
        }
        
        override func skillImpl(_ p: Point) -> SkillEffect {
            return ReaperSkillEffect(type: TYPE_REAPER_SKILL_EFFECT, x: p.x, y: p.y, radius: REAPER_SKILL_RADIUS, duration: REAPER_SKILL_DURATION, order: REAPER_SKILL_ORDER, reaper: self);
        }
    }
    
    class Destroyer: Looter {
        init( player:Player, x:Double, y:Double ) {
            super.init(type: LOOTER_DESTROYER, player: player, x: x, y: y, skillCost: DESTROYER_SKILL_COST, skillRange: DESTROYER_SKILL_RANGE, skillActive: DESTROYER_SKILL_ACTIVE)
            mass = DESTROYER_MASS;
            friction = DESTROYER_FRICTION;
            skillActive = DESTROYER_SKILL_ACTIVE;
        }
        
        override func skillImpl(_ p: Point) -> SkillEffect{
            return DestroyerSkillEffect(type: TYPE_DESTROYER_SKILL_EFFECT, x: p.x, y: p.y, radius: DESTROYER_SKILL_RADIUS, duration: DESTROYER_SKILL_DURATION,
                                        order: DESTROYER_SKILL_ORDER, destroyer: self);
        }
    }
    
    class Doof: Looter {
        init(player:Player, x:Double, y:Double) {
            super.init(type: LOOTER_DOOF, player: player, x: x, y: y, skillCost: DOOF_SKILL_COST, skillRange: DOOF_SKILL_RANGE, skillActive: DOOF_SKILL_ACTIVE)
            mass = DOOF_MASS;
            friction = DOOF_FRICTION;
            skillActive = DOOF_SKILL_ACTIVE;
        }
        
        override func skillImpl(_ p: Point) -> SkillEffect {
            return DoofSkillEffect(type: TYPE_DOOF_SKILL_EFFECT, x: p.x, y: p.y, radius: DOOF_SKILL_RADIUS, duration: DOOF_SKILL_DURATION, order: DOOF_SKILL_ORDER, doof: self);
        }
        
        // With flame effects! Yeah!
        func sing() -> Int {
            return Int(floor(speed() * DOOF_RAGE_COEF))
        }
    }
    
    
    class TankerSpawn {
        let size: Int
        let angle: Double
        
        init(size: Int, angle: Double) {
            self.size = size;
            self.angle = angle;
        }
    }
    
    class Player {
        var score = 0
        let index: Int
        var rage = 0
        var looters : [Looter] = []
        var dead : Bool = false
        var tankers: [TankerSpawn]  = []
        
        init( index: Int ) {
            self.index = index;
        }
        
        func kill() {
            dead = true;
        }
        
        func getReaper() -> Reaper? {
            guard looters.count > LOOTER_REAPER else {
                return nil
            }
            return looters[LOOTER_REAPER] as? Reaper
        }
        
        func getDestroyer() -> Destroyer? {
            guard looters.count > LOOTER_DESTROYER else {
                return nil
            }
            return looters[LOOTER_DESTROYER] as? Destroyer
        }
        
        func getDoof() -> Doof? {
            guard looters.count > LOOTER_DOOF else {
                return nil
            }
            return looters[LOOTER_DOOF] as? Doof
        }
    }
 
    class Collision {
        let t: Double
        let a: Unit?
        let b: Unit?
        
        init(t:Double, a:Unit?, b:Unit?) {
            self.t = t
            self.a = a
            self.b = b
        }
        
        convenience init(t:Double) {
            self.init(t:t, a:nil, b:nil)
        }
        
        convenience init(t:Double, a:Unit) {
            self.init(t:t, a:a, b:nil)
        }
        
        func dead() -> Tanker? {
            guard let a = a, let b = b else { return nil }
            
            if a is Destroyer && b is Tanker {
                let tanker = b as! Tanker
                if tanker.mass < REAPER_SKILL_MASS_BONUS {
                    return tanker
                }
            }
            
            if a is Tanker && b is Destroyer {
                let tanker = a as! Tanker
                if tanker.mass < REAPER_SKILL_MASS_BONUS {
                    return tanker
                }
            }
            return nil
        }
        
    }
    
    class SkillEffect : Hashable {
        
        let id: Int
        let type: Int
        var pos: Point
        let radius: Double
        var duration: Int
        var order: Int
        var known: Bool = false
        let looter: Looter
        
        init(type:Int, x:Double, y:Double, radius:Double, duration:Int, order:Int, looter:Looter) {
            id = GLOBAL_ID.next()!
            pos = Point(x,y)
            self.type = type
            self.radius = radius
            self.duration = duration
            self.looter = looter
            self.order = order
        }
        
        func apply(_ units: [Unit]) {
            duration -= 1;
            applyImpl(units.filter({ pos.isInRange($0.pos, radius + $0.radius) }))
        }
        
        func toFrameData() -> String {
            if (known) {
                return "\(id)"
            }
            
            known = true;
            return "\(id) \(round(pos.x)) \(round(pos.y)) \(looter.id) 0 \(type) \(round(radius))"
        }
        
        func applyImpl(_ units: [Unit]) {
            // abstract func -- do nothing
        }
        
        var hashValue: Int {
            let prime = 31
            let result = 1
            return prime * result + id
        }
        
        static func ==(lhs: Ref.SkillEffect, rhs: Ref.SkillEffect) -> Bool {
            if lhs === rhs {
                return true
            }
            if lhs.id != rhs.id {
                return false
            }
            return true
        }
    }

    class ReaperSkillEffect : SkillEffect {
        init(type:Int, x:Double, y:Double, radius:Double, duration:Int, order:Int, reaper:Reaper) {
            super.init(type:type, x:x, y:y, radius:radius, duration:duration, order:order, looter:reaper)
        }
        
        override func applyImpl(_ units:[Unit]) {
            // Increase mass
            for unit in units {
                unit.mass += REAPER_SKILL_MASS_BONUS
            }
        }
    }
    
    class DestroyerSkillEffect : SkillEffect {
        init(type:Int, x:Double, y:Double, radius:Double, duration:Int, order:Int, destroyer:Destroyer) {
            super.init(type:type, x:x, y:y, radius:radius, duration:duration, order:order, looter:destroyer)
        }

        override func applyImpl(_ units:[Unit]) {
            // Increase mass
            for unit in units {
                unit.thrust(pos, -DESTROYER_NITRO_GRENADE_POWER)
            }
        }
    }
    
    class DoofSkillEffect : SkillEffect {
        init(type:Int, x:Double, y:Double, radius:Double, duration:Int, order:Int, doof:Doof) {
            super.init(type:type, x:x, y:y, radius:radius, duration:duration, order:order, looter:doof)
        }
        
        override func applyImpl(_ units:[Unit]) {
            // Nothing to do now
        }
    }


    let playerCount: Int
    var units: [Unit] = []
    var looters: [Looter] = []
    var tankers: [Tanker] = []
    var deadTankers: [Tanker] = []
    var wrecks: [Wreck] = []
    var players: [Player] = []
    // var frameData: [String]
    var skillEffects: Set<SkillEffect> = []
    var random: Random

    
    func spawnTanker(_ player: Player) {
        let spawn = player.tankers.removeFirst()
        let angle = (Double(player.index) + spawn.angle) * Double.pi * 2.0 / Double(playerCount)
        let _cos = cos(angle)
        let _sin = sin(angle);
        
        if (SPAWN_WRECK) {
            // Spawn a wreck directly
            let wreck = Wreck(x:_cos * WATERTOWN_RADIUS, y:_sin * WATERTOWN_RADIUS, water:spawn.size,
                              radius:TANKER_RADIUS_BASE + Double(spawn.size) * TANKER_RADIUS_BY_SIZE);
            wreck.player = player
            wrecks.append(wreck)
            return
        }
        
        let tanker = Tanker(size: spawn.size, player: player);
        var distance = TANKER_SPAWN_RADIUS + tanker.radius;
        
        safetyCheck: while true {
            tanker.pos = Point(_cos * distance, _sin * distance)
            for unit in units {
                if tanker.pos.distance(to:unit.pos) <= tanker.radius + unit.radius {
                    distance += TANKER_MIN_RADIUS
                    continue safetyCheck
                }
            }
            break
        }
        
        tanker.thrust(WATERTOWN, TANKER_START_THRUST);
        
        units.append(tanker);
        tankers.append(tanker);
    }
    
    func createLooter(_ type: Int, _ player: Player, _ x:Double, _ y:Double) -> Looter? {
        if (type == LOOTER_REAPER) {
            return Reaper(player:player, x:x, y:y)
        } else if (type == LOOTER_DESTROYER) {
            return Destroyer(player:player, x:x, y:y)
        } else if (type == LOOTER_DOOF) {
            return Doof(player:player, x:x, y:y)
        }
        return nil
    }

    init(playerCount: Int, seed: UInt64?) throws {
        self.playerCount = playerCount;
        
        if let seed = seed {
            random = Random(seed:seed)
        } else {
            random = Random()
        }
        
        TANKERS_BY_PLAYER = TANKERS_BY_PLAYER_MIN + random.int(to:TANKERS_BY_PLAYER_MAX - TANKERS_BY_PLAYER_MIN + 1)
        
        // Create players
        for i in 0..<playerCount {
            players.append(Player(index: i))
        }
        
        // Generate the map
        var queue: [TankerSpawn] = []
        for _ in 0..<500 {
            queue.append(TankerSpawn(size: TANKER_MIN_SIZE + random.int(to:TANKER_MAX_SIZE - TANKER_MIN_SIZE), angle: random.double()))
        }
        for player in players {
            player.tankers.append(contentsOf: queue)
        }
        
        // Create looters
        for player in players {
            for i in 0..<LOOTER_COUNT {
                if let looter = createLooter(i, player, 0, 0) {
                    player.looters.append(looter)
                    units.append(looter)
                    looters.append(looter)
                }
            }
        }
        
        // Random spawns for looters
        var finished = false;
        while (!finished) {
            finished = true;
            
            for i in 0..<LOOTER_COUNT {
                if !finished {
                    break
                }
                let distance = random.double() * (MAP_RADIUS - LOOTER_RADIUS);
                let angle = random.double();
                
                for player in players {
                    let looterAngle = (Double(player.index) + angle) * (Double.pi * 2.0 / Double(playerCount));
                    let _cos = cos(looterAngle)
                    let _sin = sin(looterAngle)
                    
                    let looter = player.looters[i];
                    looter.pos = Point(_cos * distance, _sin * distance)
                    
                    // If the looter touch a looter, reset everyone and try again
                    for unit in units where unit != looter {
                        if looter.pos.distance(to: unit.pos) <= looter.radius + unit.radius {
                            finished = false
                            for l in looters {
                                l.pos = Point(0,0)
                            }
                            break
                        }
                    }
                }
            }
        }
        
        // Spawn start tankers
        for _ in 0..<TANKERS_BY_PLAYER {
            for player in players {
                spawnTanker(player)
            }
        }
        
        adjust()
        // newFrame(1.0);
        // snapshot();
    }
    
    func getInitInputForPlayer(_ playerIdx: Int) -> [String] {
        // No init input
        return []
    }
    
    func prepare(_ round:Int) {
        for looter in looters {
            looter.reset()
        }
    }

    
    func getPlayerId(_ id:Int, _ forId:Int) -> Int {
        // This method can be called with id=-1 because of the default player for units
        if (id < 0) {
            return id;
        }
        
        if (id == forId) {
            return 0;
        }
        
        if (id < forId) {
            return id + 1;
        }
        
        return id;
    }
    
    func getInputForPlayer(_ turn:Int, _ playerIdx:Int) -> [String] {
        var lines: [String] = []
        
        // Scores
        // My score is always first
        lines.append("\(players[playerIdx].score)")
        for i in 0..<playerCount where i != playerIdx {
            lines.append("\(players[i].score)")
        }
        
        // Rages
        // My rage is always first
        lines.append("\(players[playerIdx].rage)")
        for i in 0..<playerCount where i != playerIdx {
            lines.append("\(players[i].rage)")
        }
        
        // Units
        var unitsLines: [String] = []
        // Looters and tankers
        unitsLines.append(contentsOf: units.map({ (u) -> String in
            return "\(u.id) \(u.type) \(getPlayerId(u.getPlayerIndex(), playerIdx)) \(u.mass) \(Int(round(u.radius))) \(Int(round(u.pos.x))) \(Int(round(u.pos.y))) \(Int(round(u.vx))) \(Int(round(u.vy))) \(u.getExtraInput()) \(u.getExtraInput2())"
        }))
        
        // Wrecks
        unitsLines.append(contentsOf: wrecks.map({ (w) -> String in
            return "\(w.id) \(TYPE_WRECK) -1  -1  \(Int(round(w.radius))) \(Int(round(w.pos.x))) \(Int(round(w.pos.y))) 0 0  \(w.water) -1"
        }))

        // Skill effects
        unitsLines.append(contentsOf: skillEffects.map({ (s) -> String in
            return "\(s.id) \(s.type) -1 -1 \(Int(round(s.radius))) \(Int(round(s.pos.x))) \(Int(round(s.pos.y))) 0 0 \(s.duration) -1"
        }))
        
        lines.append("\(unitsLines.count)")
        lines.append(contentsOf: unitsLines)
        
        return lines
    }
    
    func getExpectedOutputLineCountForPlayer(_ playerIdx:Int) -> Int {
        return 3
    }
    
    enum Action {
        case SKILL, MOVE, WAIT
    }
    
    enum SkillResultCode {
        case OK, NO_RAGE, TOO_FAR
    }
    
    struct SkillResult {
        let target: Point
        var code: SkillResultCode
        
        init(_ target: Point) {
            self.target = target
            code = .OK
        }
        
        func getX() -> Int {
            return Int(target.x)
        }
        
        func getY() -> Int {
            return Int(target.y)
        }
    }


    enum RefException: Error {
        case WinException, LostException, InvalidInputException, GameOverException
    }
    
    func handlePlayerOutput(_ frame:Int, _ turn: Int, _ playerIdx: Int, _ outputs: [Move]) {
        
        //  BUGBUG -- for simplicity (and since we aren't running in containers) remove the text response and
        // replace with codes.  Also messages are getting the short end of the stick
        let player = players[playerIdx]
        for i in 0..<LOOTER_COUNT {
            let move = outputs[i]
            let looter = player.looters[i]
            switch move {
            case .wait:
                looter.attempt = .WAIT
                continue
            case .accel(let target, let power):
                looter.setWantedThrust(target, power)
                looter.attempt = .MOVE
                continue
            case .skill(let target):
                if !looter.skillActive {
                    // Don't kill the player for that. Just do a WAIT instead
                    looter.attempt = Action.WAIT
                    continue
                }

                looter.attempt = .SKILL
                var result = SkillResult(target)
                
                do {
                    let effect = try looter.skill(target)
                    skillEffects.insert(effect!);
                } catch LooterError.NoRageException {
                    result.code = .NO_RAGE
                } catch LooterError.TooFarException {
                    result.code = .TOO_FAR
                } catch {
                    // generic catchall -- no other exceptions
                }
                looter.skillResult = result
                continue
            }
        }
    }
    

    // Get the next collision for the current round
    // All units are tested
    func getNextCollision() -> Collision {
        var result:Collision = Ref.NULL_COLLISION
        
        for i in 0..<units.count {
            let unit = units[i]
            // Test collision with map border first
            let collision = unit.getCollision()
            if collision.t < result.t {
                result = collision
            }
            
            for j in (i+1)..<units.count {
                let collision = unit.getCollision(units[j])
                if collision.t < result.t {
                    result = collision;
                }
            }
        }
        
        return result;
    }
    
    // Play a collision
    func playCollision(_ collision:Collision) {
        guard let a = collision.a else { return }
        if let b = collision.b {
            if let dead = collision.dead() {
                // A destroyer kill a tanker
                // addDeadToFrame(dead);
                deadTankers.append(dead);
                if let index = tankers.index(of: dead) {
                    tankers.remove(at: index)
                }
                if let index = units.index(of: dead) {
                    units.remove(at: index)
                }
                
                if let wreck = dead.die() {
                    // If a tanker is too far away, there's no wreck
                    wrecks.append(wreck);
                    // addToFrame(wreck);
                }
            } else {
                // Bounce between two units
                //addToFrame(collision.a);
                //addToFrame(collision.b);
                a.bounce(b)
            }
        } else {
            // Bounce with border
            // addToFrame(collision.a);
            a.bounce();
        }
    }
    
    
    func updateGame(_ turn:Int) {
        // Apply skill effects
        for effect in skillEffects {
            effect.apply(units)
        }
    
        // Apply thrust for tankers
        for tanker in tankers {
            tanker.play();
        }
        
        // Apply wanted thrust for looters
        for player in players {
            for looter in player.looters {
                if let target = looter.wantedThrustTarget {
                    looter.thrust(target, looter.wantedThrustPower)
                }
            }
        }
        
        var t = 0.0;
    
        // Play the round. Stop at each collisions and play it. Repeat until t > 1.0
        var collision = getNextCollision();
        while (collision.t + t <= 1.0) {
            let delta = collision.t;
            for unit in units {
                unit.move(delta)
            }
            t += collision.t;
            playCollision(collision);
            collision = getNextCollision();
        }
        
        // No more collision. Move units until the end of the round
        let delta = max(0,1.0 - t)
        for unit in units {
            unit.move(delta)
        }

        var tankersToRemove:[Tanker] = []

        for tanker in tankers {
            let distance = tanker.pos.distance(to:WATERTOWN)
            let full = tanker.isFull();
            
            if (distance <= WATERTOWN_RADIUS && !full) {
                // A non full tanker in watertown collect some water
                tanker.water += 1;
                tanker.mass += TANKER_MASS_BY_WATER;
            } else if (distance >= TANKER_SPAWN_RADIUS + tanker.radius && full) {
                // Remove too far away and not full tankers from the game
                tankersToRemove.append(tanker);
            }
        }
    
        for tanker in tankersToRemove {
            if let index = units.index(of: tanker) {
                units.remove(at: index)
            }
            if let index = tankers.index(of: tanker) {
                tankers.remove(at: index)
            }
            deadTankers.append(tanker)
        }

        // Spawn new tankers for each dead tanker during the round - including
        // those killed by collision above
        for tanker in deadTankers {
            spawnTanker(tanker.player)
        }
        deadTankers = []
        
        var deadWrecks = Set<Wreck>()
        
        wrecks = wrecks.filter( {(w) in
            let alive = w.harvest(players: players, skillEffects: skillEffects)
            if !alive {
                //addDeadToFrame(w)
                deadWrecks.insert(w)
            }
            return alive
        })
        
        if (SPAWN_WRECK) {
            for dw in deadWrecks {
                spawnTanker(dw.player!)
            }
        }
        
        // Round values and apply friction
        adjust();
        
        // Generate rage
        if (LOOTER_COUNT >= 3) {
            for player in players {
                if let doof = player.getDoof() {
                    player.rage = min(MAX_RAGE, player.rage + doof.sing())
                }
            }
        }
        
        // Restore masses
        for unit in units {
            while unit.mass >= REAPER_SKILL_MASS_BONUS {
                unit.mass -= REAPER_SKILL_MASS_BONUS
            }
        }
        
        // Remove dead skill effects
        var effectsToRemove:Set<SkillEffect> = []
        for se in skillEffects {
            if se.duration <= 0 {
                // addDeadToFrame(ef)
                effectsToRemove.insert(se)
            }
        }
        for se in effectsToRemove {
            skillEffects.remove(se)
        }
    }
    
    func adjust() {
        for unit in units {
            unit.adjust(skillEffects: skillEffects)
        }
    }
    
    func isGameOver() -> Bool {
        for player in players {
            if player.score >= WIN_SCORE {
                return true
            }
        }

        let alive = players.filter {!$0.dead}
        if alive.count == 1 {
            let survivor = alive[0]
            for player in players where player.index != survivor.index {
                if player.score > survivor.score {
                    return false
                }
            }
            return true
        }
        return alive.isEmpty
    }
    
    func getWinner() -> Int? {
        var highScore = 0
        var winner = -1
        for player in players {
            if player.score > highScore && player.dead == false {
                highScore = player.score
                winner = player.index
            }
        }
        if winner == -1 {
            return nil
        }
        return winner
    }
    
    func getMinimumPlayerCount() -> Int {
        return 3
    }
    
    func getScore(_ playerIdx: Int) -> Int {
        return players[playerIdx].score
    }

    func getMaxRoundCount(_ playerCount: Int) -> Int {
        return 200
    }
    
    func getMsTime(for round:Int) -> Int {
        if round == 1 {
            return 1000
        }
        return 50
    }
    
    func displayState(for turn:Int) {
        
        let xScale = 500
        let yScale = 600
        
        debugPrint( "****************************************************************************************************", to:&errStream )
        debugPrint( " Round: \(turn) -- P1: \(players[0].score) -- P2: \(players[1].score) -- P3: \(players[2].score)", to: &errStream)
        debugPrint( "****************************************************************************************************", to:&errStream )
        let minY = Int(MAP_RADIUS)/yScale * -1
        let maxY = Int(MAP_RADIUS)/yScale
        let minX = Int(MAP_RADIUS)/xScale * -1
        let maxX = Int(MAP_RADIUS)/xScale
        for y in minY...maxY {
            var r1 = ""
            var r2 = ""
            var r3 = ""
            for x in minX...maxX {
                
                var s1: String = ""
                var s2: String = ""
                var s3: String = ""
                
                for unit in units {
                    if Int(round(unit.pos.x))/xScale == x &&  Int(round(unit.pos.y))/yScale == y {
                        s1 += UnitType(unit.type).shortName
                        s2 += String(unit.getPlayerIndex()+1)
                    }
                }
                while s1.count < 4 {
                    s1 += " "
                    s2 += " "
                }

                for w in wrecks {
                    if Int(round(w.pos.x))/xScale == x &&  Int(round(w.pos.y))/yScale == y {
                        s3 += UnitType.wreck.shortName
                    }
                }
                for s in skillEffects {
                    if Int(round(s.pos.x))/xScale == x &&  Int(round(s.pos.y))/yScale == y {
                        s3 += UnitType(s.type).shortName
                    }
                }
                while s3.count < 4 {
                    s3 += " "
                }

                r1 += s1
                r2 += s2
                r3 += s3
            }
            debugPrint( r1, to:&errStream )
            debugPrint( r2, to:&errStream )
            debugPrint( r3, to:&errStream )
        }
    }

}

protocol PlayerAgent {
    func initialize()
    func initialize( withData data: [String] )
    func takeTurn() -> [Move]
    func takeTurn( withData data: [String] ) -> [Move]
}

class GenericAgent : PlayerAgent {
    init() {
    }
    
    func initialize() {
        parseInitialInput(rawData: readInitialInput())
    }
    
    func initialize(withData data: [String]) {
        parseInitialInput(rawData: data)
    }
    
    func takeTurn() -> [Move] {
        return takeTurn( withData: readInput() )
    }
    
    func takeTurn( withData data: [String] ) -> [Move] {
        return [Move.wait, Move.wait, Move.wait]
    }
    
    func readInitialInput() -> [String] {
        return [] // no initial input
    }
    
    func parseInitialInput( rawData data: [String] ) {
        return // no initial input
    }
    
    func readInput() -> [String] {
        var data: [String] = []
        
        // Scores
        data.append(readLine()!)
        data.append(readLine()!)
        data.append(readLine()!)
        
        // Rage(s)
        data.append(readLine()!)
        data.append(readLine()!)
        data.append(readLine()!)
        
        let line = readLine()!
        data.append(line)
        let unitCount = Int(line)!
        if unitCount > 0 {
            for _ in 0..<unitCount {
                data.append(readLine()!)
            }
        }
        
        return data
    }
    
    func parseInput( rawData data: [String] ) {
        var iter = data.makeIterator()
        
        let p0Score = Int(iter.next()!)!
        let p1Score = Int(iter.next()!)!
        let p2Score = Int(iter.next()!)!
        
        let p0Rage = Int(iter.next()!)!
        let p1Rage = Int(iter.next()!)!
        let p2Rage = Int(iter.next()!)!
        
        let unitCount = Int(iter.next()!)!
        while let line = iter.next() {
            let inputs = line.characters.split{$0 == " "}.map(String.init)
            let unitId = Int(inputs[0])!
            let unitType = UnitType( Int(inputs[1])! )
            let player = Int(inputs[2])!
            let mass = Double(inputs[3])!
            let radius = Double(inputs[4])!
            let x = Int(inputs[5])!
            let y = Int(inputs[6])!
            let vx = Int(inputs[7])!
            let vy = Int(inputs[8])!
            let extra = Int(inputs[9])!
            let extra2 = Int(inputs[10])!
        }
    }
}

class Game {
    
    var ref: Ref
    let playerCount = 3
    
    init( seed: UInt64? = nil ) {
        ref = try! Ref(playerCount: playerCount, seed: seed)
    }
    
    func run() {
        
        var turn = 0
        var timer = Timer()

        var agents = [ GenericAgent(), GenericAgent(), GenericAgent() ]
        
        while turn < ref.getMaxRoundCount(playerCount) {
            
            turn += 1
            
            ref.displayState(for: turn)
            ref.prepare(turn)
            
            var playerInput: [[String]] = []
            for i in 0..<playerCount {
                let data = ref.getInputForPlayer(turn, i)
                playerInput.append(data)
            }

            var playerMoves: [[Move]] = []
            for i in 0..<playerCount {
                timer.start()
                let moves = agents[i].takeTurn(withData: playerInput[i])
                let t = timer.peek_mSec()
                if t > Double(ref.getMsTime(for:turn)) {
                    // technically the player should lose here for taking too much time, but since we're testing
                    // we'll spit out a warning instead
                    debugPrint("Player \(i) is over time at \(t) mSec")
                }
                playerMoves.append(moves)
            }
            
            for i in 0..<playerCount {
                debugPrint("\(i): {\(playerMoves[i].map({$0.description}).reduce("", {"\($0) \($1)"})) }", to: &errStream)
                let _ = ref.handlePlayerOutput(0, turn, i, playerMoves[i])
            }
            
            ref.updateGame(turn)
            
            if ref.isGameOver() {
                break
            }
        }
        
        if let winner = ref.getWinner() {
            debugPrint("\(winner + 1) wins!", to: &errStream)
        } else {
            debugPrint("No winners!", to: &errStream)
        }
    }
}

