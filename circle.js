function sphere(x, y, z, radius){
    this.x = x;
    this.y = y;
    this.z = z;
    this.radius = radius;
}

function rand(min, max) {
    return Math.random() * (max - min) + min;
}

function sphereToSphereCollision(sphere1, sphere2) {
    return (sphere1.radius + sphere2.radius) > Math.hypot(Math.hypot(sphere1.x - sphere2.x, sphere1.y - sphere2.y), sphere1.z - sphere2.z)
}

function reassignSphere(sphere){
    sphere.x = rand(-1000, 1000);
    sphere.y = rand(-1000, 1000);
    sphere.z = rand(-1000, 1000);
    sphere.radius = rand(-1000, 1000);
}

function runSphereToSphereCollisionJS(loops) {
    var spheres = []
    for (var i = 0; i < loops; i++)
    {
        spheres.push(new sphere(rand(-1000, 1000), rand(-1000, 1000), rand(-1000, 1000), rand(-1000, 1000)))
    }
    var t0 = performance.now()
    var total = 0.0;
    for (var i = 0; i < loops; i+=2)
    {
        total += sphereToSphereCollision(spheres[i], spheres[i+1]);
    }
    var t1 = performance.now()
    console.log('total was ' + total);
    console.log('duration was ' + (t1 - t0));
    document.getElementById('wasm-js').innerText = 'duration was ' + (t1 - t0);
}
