{
  description = "Development shell for babyjubjub-ec";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs =
    { nixpkgs, ... }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];

      forAllSystems = nixpkgs.lib.genAttrs systems;
    in
    {
      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          rustTools = [
            pkgs.rustc
            pkgs.cargo
            pkgs.clippy
            pkgs.cargo-audit
          ];
        in
        {
          default = pkgs.mkShell {
            packages = rustTools;
          };
        }
      );

      apps = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          runtimeInputs = [
            pkgs.rustc
            pkgs.cargo
            pkgs.clippy
            pkgs.cargo-audit
          ];
          clippy = pkgs.writeShellApplication {
            name = "babyjubjub-ec-clippy";
            inherit runtimeInputs;
            text = ''
              if [ "$#" -eq 0 ]; then
                set -- --locked --all-features --all-targets -- -D warnings
              fi

              exec cargo clippy "$@"
            '';
          };
          audit = pkgs.writeShellApplication {
            name = "babyjubjub-ec-audit";
            inherit runtimeInputs;
            text = ''
              if [ "$#" -eq 0 ]; then
                set -- audit
              fi

              exec cargo-audit "$@"
            '';
          };
        in
        {
          clippy = {
            type = "app";
            program = "${clippy}/bin/babyjubjub-ec-clippy";
            meta.description = "Run cargo clippy for babyjubjub-ec";
          };

          audit = {
            type = "app";
            program = "${audit}/bin/babyjubjub-ec-audit";
            meta.description = "Run cargo audit for babyjubjub-ec";
          };
        }
      );
    };
}
