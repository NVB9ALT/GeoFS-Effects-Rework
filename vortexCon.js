function runVortexCons() {
geofs.fx.Particle.prototype = {
    create: function () {
        this._currentScale = this._options.startScale;
        this._currentOpacity = this._options.startOpacity;
        this._currentRotation = this._options.startRotation;
        this._options = Object.assign(this._options, { opacity: this._currentOpacity, scale: this._currentScale, rotation: this._currentRotation, color: geofs.fx.cloudManager.cloudColor }, geofs.fx.particleBillboardOptions);
        this._options.model
            ? ((this._options.rotation = [this._currentRotation, 0, 0]), (this._APIElement = new geofs.api.Model(this._options.model, this._options)))
            : this._options.groundTexture
            ? (this._APIElement = new geofs.api.groundTexture(this.currentLocation, this._options.groundTexture, this._options))
            : ((this._options.distanceDisplayCondition = new Cesium.DistanceDisplayCondition(this._options.near || 1, this._options.far || 1e6)),
              (this._APIElement = new geofs.api.billboard(this.currentLocation, this._options.url, this._options)));
    },
    setColor: function (a) {
        this._APIElement && this._APIElement.setColor(a);
    },
    setLocation: function (a) {
        this._APIElement && ((this.currentLocation = a), this._APIElement.setLocation(this.currentLocation));
    },
    setRotation: function (a, b) {
        this._APIElement && ((this._currentRotation = a), this._APIElement.setRotation(this._currentRotation, b));
    },
    setScale: function (a) {
        this._APIElement && this._APIElement.setScale(a);
    },
    setPositionOrientationAndScale: function (a, b, c) {
        this._APIElement && this._APIElement.setPositionOrientationAndScale(a, b, c);
    },
    update: function (a) {
        var b = (geofs.utils.fastNow() - this._birth) / this._options.life;
        1 < b
            ? this.destroy()
            : this._APIElement &&
              ((b = geofs.utils.easingFunctions[this._options.easing](b)),
              this._options.dtOpacity && ((this._currentOpacity = clamp(this._options.startOpacity + this._options.dtOpacity * b, 0, 1)), this._APIElement.setOpacity(this._currentOpacity)),
              this._options.dtScale && ((this._currentScale = this._options.startScale + this._options.dtScale * b), this._APIElement.setScale(this._currentScale)),
              this._options.dtRotation && ((this._currentRotation = this._options.startRotation + this._options.dtRotation * b), this._APIElement.setRotation(this._currentRotation)),
              this._options.velocity &&
                  this._options.direction &&
                  (this._options.velocityDamper && (this._options.velocity = this._options.velocity * this._options.velocityDamper * a),
                  (this.currentLocation = V3.add(this.currentLocation, V3.scale(this._options.direction, this._options.velocity))),
                  this._APIElement.setLocation(this.currentLocation)),
              this._options.floatsOnWaves && ((this.currentLocation[2] = geofs.fx.water.getWaveHeight(this.currentLocation[0], this.currentLocation[1])), this._APIElement.setLocation(this.currentLocation)));
    },
    destroy: function () {
        this._APIElement && (this._APIElement.destroy(), (this._APIElement = null));
        this._emitter = null;
        geofs.fx.particles[this._id] = null;
        delete geofs.fx.particles[this._id];
    },
};

//not using regular variables because AFAIK you can't use them to set the worldPosition of a particle emitter
geofs.bodyID = 0;
geofs.tipIndexL = 0;
geofs.tipIndexR = 0;
geofs.VCONarray = [];
var condensing = new Boolean(0);

updateSourcePoints = function() {
geofs.aircraft.instance.definition.parts.forEach(function(a, b){
   if (a.name.includes("body") || a.name.includes("Body")) {
geofs.bodyID = b
geofs.VCONarray = a.collisionPoints
geofs.VCONarray.forEach(function(c, d){
   if (c[0] > geofs.VCONarray[geofs.tipIndexR][0]){ geofs.tipIndexR = d; };
   if (c[0] < geofs.VCONarray[geofs.tipIndexL][0]){ geofs.tipIndexL = d; };
})
	}
})
}

function updateVCondensation() {
   if (geofs.cons == true && (geofs.animation.values.aoa > 10 || geofs.animation.values.loadFactor > 4) && geofs.aircraft.instance.definition.mass > 10000 && condensing != 1) {
geofs.fx.vcondensationEmitterLeft = new geofs.fx.ParticleEmitter({
				anchor: geofs.aircraft.instance.definition.parts[geofs.bodyID].collisionPoints[geofs.tipIndexL],
            duration: 1E10,
            rate: .5,
            life: 1000,
            easing: "easeOutQuart",
            startScale: .0015,
            endScale: .0015,
            randomizeStartScale: 0.00001,
            randomizeEndScale: .0001,
            startOpacity: 0.9,
            endOpacity: 0.0001,
            startRotation: "random",
            texture: "whitesmoke"
})
geofs.fx.vcondensationEmitterRight = new geofs.fx.ParticleEmitter({
				anchor: geofs.aircraft.instance.definition.parts[geofs.bodyID].collisionPoints[geofs.tipIndexR],
            duration: 1E10,
            rate: .5,
            life: 1000,
            easing: "easeOutQuart",
            startScale: .0015,
            endScale: .0015,
            randomizeStartScale: 0.00001,
            randomizeEndScale: .0001,
            startOpacity: 0.9,
            endOpacity: 0.0001,
            startRotation: "random",
            texture: "whitesmoke"
})
condensing = 1
   } else if (condensing == 1 && (geofs.animation.values.aoa < 10 || geofs.animation.values.loadFactor < 4 || geofs.aircraft.instance.definition.mass < 10000)) {
geofs.fx.vcondensationEmitterLeft.destroy()
geofs.fx.vcondensationEmitterRight.destroy()
condensing = 0
	}
}
condensationInt = setInterval(function(){
   updateSourcePoints()
	updateVCondensation()
},100)
}
